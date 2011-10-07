(ns clj-mcmc.core
  "clj-mcmc.core provides the core functions for sampling with MCMC algorithms."
  )

(defn- test-candidate [old cand propose-dens true-dens u]
  (let [denom (* (true-dens old) (propose-dens cand))]
    (if (and (not= denom 0.)
             (> (/ (* (true-dens cand) (propose-dens old)) denom) u))  
      cand
      old)))

(defn mcmc [start propose-fn propose-dens true-dens rn-stream]
  (cons start
        (lazy-seq
         (let [rn1 (first rn-stream)
               rn2 (second rn-stream)
               rn-next (nthnext rn-stream 2)]
           (mcmc
            (test-candidate start (propose-fn start rn1) 
                            propose-dens true-dens (to-uniform rn2))
            propose-fn propose-dens true-dens rn-next)))))
