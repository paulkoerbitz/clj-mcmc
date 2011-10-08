(ns clj-mcmc.distributions
  "Alternative functions for those that we don't reuse from Incanter
   (mostly sampling functions)"
  (:use (clj-mcmc helpers))
  (:require [incanter.core :as ic]
            [incanter.stats :as is]))

(defn to-uniform [rand & {:keys [min max] :or {min 0.0 max 1.0}}]
  (+ (* rand (- max min)) min))

(defn to-normal [rand & {:keys [mu std] :or {mu 0.0 std 1.0}}]
  (+ (* (is/quantile-normal rand) std) mu))

(defn to-multivar-normal [rnd &
                          {:keys [mu cov]
                           :or {mu (ic/matrix (replicate (count rnd) 0)) 
                                cov (ic/identity-matrix (count rnd))}}]
  (ic/plus (ic/mmult cov (ic/matrix (is/quantile-normal rnd))) mu))


(defn pdf-multivar-normal [x & {:keys [mu cov]
                                :or {mu (ic/matrix (replicate (count x) 0))
                                     cov (ic/identity-matrix (count x))}}]
  (let [n (count x)
        inv-cov (ic/solve cov)
        x-mu (ic/minus (ic/matrix x) mu)]
    (/ (Math/exp (/ (ic/mmult (ic/mmult (ic/trans x-mu) inv-cov) x-mu)  -2.0))
       (* (Math/pow (* 2.0 Math/PI) (/ n 2.0))
          (Math/sqrt (ic/det cov))))))
