(ns examples.asset-model.log-likelihood)

(defn black-scholes-vasicek-loglik
  "Loglikelihood function of the black-scholes-vasicek model"
  [[mu sS ka be sR ro] ^floats stocks ^floats rs dt]
  (let [n (int (alength stocks))
        ekdt (Math/exp (- (* ka dt)))
        t_1_ekdt (/ (* be (- 1.0 ekdt)) ka)
        sig_ls  (* sS (Math/sqrt dt))
        sig_ls_sq (* sig_ls sig_ls)
        sig_r (* sR (Math/sqrt (/ (- 1.0 (* ekdt ekdt)) (* 2.0 ka))))  
        sig_r_sq (* sig_r sig_r)
        rho_sq (* ro ro)
        drift_ls (- (* mu dt) (/ sig_ls_sq 2.0))
        inner-sum
        (loop [i (int 1), sum (float 0.0),
               log_s_before (float (Math/log (aget stocks 0)))]
          (if (< i (alength stocks))
            (let [log_s_now (float (Math/log (aget stocks i))) 
                  x1 (float (- log_s_now log_s_before drift_ls)) 
                  x2 (float (- (aget rs i) (* (aget rs (dec i)) ekdt) t_1_ekdt))]
              (recur (inc i)
                     (float (+ sum
                               (/ (* x1 x1) sig_ls_sq)
                               (- (/ (* 2.0 ro x1 x2) (* sig_ls sig_r)) )
                               (/ (* x1 x1) sig_r_sq)))
                     log_s_now))
            sum))]
    (/ (+ (* n (Math/log (* (- 1.0 rho_sq) sig_ls_sq sig_r_sq)))
          (/ inner-sum (- 1.0 rho_sq))) 
       -2.0)))

(defn cev-ckls-loglik
  "log-likelihood function for the cev-ckls model"
  [[mu sS al ka be sR xi ro] ^floats stocks ^floats rs dt]
  (if (or (< sS 1e-10) (< al 1e-10) (< sR 1e-10) (< xi 1e-10))
    -1e308
    (let [rho_comp (float (- 1.0 (* ro ro)))
          sqrtDt (float (Math/sqrt dt))]
      (loop [i (int 1), sum (float 0.0)]
        (if (< i (alength stocks))
          (let [x (float (- (aget stocks i)
                            (* (+ 1.0 (* mu dt)) (aget stocks (dec i)))))   
                y (float (- (aget rs i) (* (aget rs (dec i)) (- 1.0 (* ka dt))) (* be dt)))
                eta_S (float (* sS (Math/pow (aget stocks (dec i)) al) sqrtDt))  
                eta_r (float (* sR (Math/pow (aget rs (dec i)) xi) sqrtDt))]
            (recur (inc i)
                   (float (- sum
                             (Math/log (* 2.0 Math/PI (Math/sqrt rho_comp) eta_S eta_r))
                             (/ (+ (/ (* x x) (* eta_S eta_S))
                                   (/ (* -2.0 ro x y) (* eta_S eta_r))
                                   (/ (* y y) (* eta_r eta_r)))
                                (* 2.0 rho_comp))))))
          sum)))))


(defn diffusion-loglik
  [f-drift f-vol ^floats xs]
  (let [n (int (count xs))]
    (loop [i     (int 1)
           x-old (aget xs 0)
           sum   (float (* 0.5 (dec n) (Math/log (* 2.0 Math/PI))))]
      (if (< i n)
        (let [x-new (aget xs i)
              x-mu  (- x-new (f-drift x-old))
              vol   (f-vol x-old)]
          (recur (inc i) x-new
                 (float (- sum (Math/log vol) (/ (* x-mu x-mu) (* 2.0 vol vol))))))
        sum))))

(defn ckls-loglik ^Double
  [[^Double ka ^Double be ^Double sR ^Double xi] ^floats rs dt]
  (let [f-drift  (fn ^Double [^Double x-old] (+ (* x-old (- 1.0 (* ka dt))) (* be dt)))
        sqrtDt  (Math/sqrt dt)
        f-vol   (fn ^Double [^Double x-old] (* sR (Math/pow x-old xi) sqrtDt))]
    (diffusion-loglik f-drift f-vol rs)))

(defn cev-loglik
  [[mu sS al] ^floats xs dt]
  (if (or (< sS 1e-6) (< al 1e-6))
    -1e308
    (let [one+mu*dt (+ 1.0 (* mu dt)) 
          f-drift   (fn [x-old] (* one+mu*dt x-old) )
          sqrtDt    (Math/sqrt dt)
          f-vol     (fn [x-old] (* sS (Math/pow x-old al) sqrtDt))]
      (diffusion-loglik f-drift f-vol xs))))
