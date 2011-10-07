(ns examples.asset-models)

(set! *warn-on-reflection* true)

(defn sq [x] (* x x))

(defn black-scholes-vasicek-loglik
  "Loglikelihood function of the black-scholes-vasicek model"
  [mu sS ka be sR ro ^floats stocks ^floats rs dt]
  (let [n (int (alength stocks))
        ekdt (Math/exp (- (* ka dt)))
        t_1_ekdt (/ (* be (- 1.0 ekdt)) ka)
        sig_ls  (* sS (Math/sqrt dt))
        sig_ls_sq (sq sig_ls)
        sig_r (/ (* sR (Math/sqrt (- 1.0 (sq ekdt)))) (* 2.0 ka)) 
        sig_r_sq (sq sig_r)
        rho_sq (sq ro)
        drift_ls (- (* mu dt) (/ sig_ls_sq 2.0))
        inner-sum
        (loop [i (int 1), sum (float 0.0),
               log_s_before (float (Math/log (aget stocks 0)))]
          (if (< i (alength stocks))
            (let [log_s_now (float (Math/log (aget stocks i))) 
                  x1 (float (- log_s_now log_s_before drift_ls)) 
                  x2 (float (aget rs i))]
              (recur (inc i)
                     (float (+ sum
                               (/ (sq x1) sig_ls_sq)
                               (- (/ (* 2.0 ro x1 x2) (* sig_ls sig_r)) )
                               (/ (sq x2) sig_r_sq)))
                     log_s_now))
            sum))]
    (/ (+ (* n (Math/log (* (- 1.0 rho_sq) sig_ls_sq sig_r_sq)))
          (/ inner-sum (- 1.0 rho_sq))) 
       -2.0)))

(defn cev-ckls-loglik
  "log-likelihood function for the cev-ckls model"
  [mu sS al ka be sR xi ro ^floats stocks ^floats rs dt]
  (let [rho_comp (float (- 1.0 (sq ro)))
        sqrtDt (float (Math/sqrt dt))]
    (loop [i (int 1), sum (float 0.0)]
      (if (< i (alength stocks))
        (let [x (float (- (aget stocks i)
                          (* (+ 1.0 (* mu dt)) (aget stocks (dec i)))))   
              y (float (- (aget rs i) (aget rs (dec i))
                          (- (* ka (aget rs (dec i)) dt) be)))
              eta_S (float (* sS (Math/pow (aget stocks (dec i)) al) sqrtDt))  
              eta_r (float (* sR (Math/pow (aget rs (dec i)) xi) sqrtDt))  
              eta_S_sq (float (sq eta_S)) 
              eta_r_sq (float (sq eta_r))]
          (recur (inc i)
                 (float
                  (- sum
                     (Math/log (* 2.0 Math/PI (sqrt rho_comp) eta_S eta_r))
                     (/ (+ (/ (sq x) eta_S_sq)
                           (/ (* -2.0 ro x y) (* eta_S eta_r))
                           (/ (sq y) eta_r_sq))
                        (* 2.0 rho_comp))))))
        sum))))

