(ns examples.simple-mcmc
"
  Showcases simple examples of how to use the MCMC algorithms
"
  (:use (clj-mcmc core distributions random helpers ))
  (:use (incanter core stats charts)))

(defn sample-from-normal
  "Samples from a normal distribution via MCMC. Samples are proposed
   from a uniform distribution and accepted or rejected based on the
   normal density."
  [n]
  (let [min-max {:min -4. :max 4.}
        samples (take n
                      (mcmc 0.0
                            (fn [prev rand] (to-uniform rand :min -4. :max 4.)) 
                            #(pdf-uniform % :min -4. :max 4.)
                            pdf-normal
                            (rand-stream)))]
    (view (histogram samples :nbins (between (int (/ n 100)) 10 1000)))))
