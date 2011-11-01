(ns clj-mcmc.core
  "clj-mcmc.core provides the core functions for sampling with MCMC algorithms."
  (:use [clj-mcmc distributions helpers]))

(defn- test-candidate [old cand propose-dens true-dens u]
  (let [denom (* (true-dens old) (propose-dens old cand))]
    (if (and (not= denom 0.)
             (> (/ (* (true-dens cand) (propose-dens old old)) denom) u))  
      cand
      old)))

(defn- test-w-logs [old cand propose-dens true-dens u]
  (if (> (- (+ (true-dens cand) (propose-dens old old))
            (+ (true-dens old) (propose-dens old cand)))
         (Math/log u))
    cand
    old))

(defn mcmc [start propose-fn propose-dens true-dens rnd-stream
            & {:keys [logs] :or {logs true}}]
  (let [n       (count start)
        test-fn (if logs test-w-logs test-candidate)]
    (cons start
        (lazy-seq
         (let [u        (to-uniform (first rnd-stream))
               rnd      (vec (take n (next rnd-stream))) 
               rnd-next (nthnext rnd-stream (+ n 1))
               cand     (propose-fn start rnd)
               next     (test-fn start cand propose-dens true-dens u)]
           
           (mcmc next propose-fn propose-dens true-dens rnd-next)))))) 

(defn thin [n coll]
  (when (seq coll)
    (cons (first coll) (lazy-seq (thin n (nthnext coll (inc n)))))))

