(ns clj-mcmc.distributions
  "Alternative functions for those that we don't reuse from Incanter
   (mostly sampling functions)"
  (:use (clj-mcmc helpers)))

(defn to-uniform [rand & {:keys [min max] :or {min 0.0 max 1.0}}]
  (+ (* rand (- max min)) min))
