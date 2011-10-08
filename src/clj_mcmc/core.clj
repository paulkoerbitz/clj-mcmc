(ns clj-mcmc.core
  "clj-mcmc.core provides the core functions for sampling with MCMC algorithms."
  (:use [clj-mcmc distributions helpers]))

(defn- test-candidate [old cand propose-dens true-dens u]
  (let [denom (* (true-dens old) (propose-dens cand))]
    (if (and (not= denom 0.)
             (> (/ (* (true-dens cand) (propose-dens old)) denom) u))  
      cand
      old)))

(defn mcmc [start propose-fn propose-dens true-dens rnd-stream]
  (let [n (count start)]
    (cons start
        (lazy-seq
         (let [u (to-uniform (first rnd-stream))
               rnd (take n (next rnd-stream))
               rnd-next (nthnext rnd-stream (+ n 1))
               cand (propose-fn start rnd)
               next (test-candidate start cand propose-dens true-dens u)
               ]
           (mcmc next propose-fn propose-dens true-dens rnd-next)))))) 
