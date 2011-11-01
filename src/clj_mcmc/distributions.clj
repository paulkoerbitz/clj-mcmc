(ns clj-mcmc.distributions
  "Alternative functions for those that we don't reuse from Incanter
   (mostly sampling functions)"
  (:use (clj-mcmc helpers))
  (:require [incanter.core :as ic]
            [incanter.stats :as is]))

(defn to-uniform [rand & {:keys [min max] :or {min 0.0 max 1.0}}]
  (+ (* rand (- max min)) min))

(defn to-normal [rand & {:keys [mean sd] :or {mean 0.0 sd 1.0}}]
  (+ (* (is/quantile-normal rand) sd) mean))

(defn to-multivar-normal [rnd &
                          {:keys [mean sd]
                           :or {mean (ic/matrix (replicate (count rnd) 0)) 
                                sd (ic/identity-matrix (count rnd))}}]
  (ic/plus (ic/mmult sd (ic/matrix (is/quantile-normal rnd))) mean))


(defn pdf-multivar-normal [x & {:keys [mu cov]
                                :or {mu (ic/matrix (replicate (count x) 0))
                                     cov (ic/identity-matrix (count x))}}]
  (let [n (count x)
        x-mu (ic/minus (ic/matrix x) mu)
        rhs (ic/solve cov x-mu)]
    (/ (Math/exp (/ (ic/mmult (ic/trans x-mu) rhs)  -2.0))
       (* (Math/pow (* 2.0 Math/PI) (/ n 2.0))
          (Math/sqrt (ic/det cov))))))

(defn to-truncated-normal [rnd & {:keys [mean sd lower-bound upper-bound]
                                  :or {mean 0.0, sd 1.0}}]
  (let [lb (if lower-bound (is/cdf-normal lower-bound) 0.0)
        ub (if upper-bound (is/cdf-normal upper-bound) 1.0)]
    (to-normal (to-uniform rnd :min lb :max ub) :mean mean :sd sd)))

(defn pdf-truncated-normal [x & {:keys [mean sd lower-bound upper-bound]
                                 :or {mean 0.0, sd 1.0}}]
  (if (or (and lower-bound (<= x lower-bound))
          (and upper-bound (>= x upper-bound)))
    0.0
    (let [l-tail (if lower-bound (is/cdf-normal lower-bound :mean mean :sd sd) 0.0)
          r-tail (if upper-bound (is/cdf-normal upper-bound :mean mean :sd sd) 1.0)]
      (/ (is/pdf-normal x :mean mean :sd sd) (- r-tail l-tail)))))
