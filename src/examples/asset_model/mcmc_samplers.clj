(ns examples.asset-model.mcmc-samplers
  (:use [clj-mcmc.helpers :only [! average]]
        [clj-mcmc.core :only [mcmc]]
        [clj-mcmc.random :only [rand-stream]])
  (:require [clj-mcmc.distributions :as dist]
            [examples.asset-model.log-likelihood :as ll]
            [incanter.core :as ic]
            [incanter.stats :as is]))

(defn sample-bsv-parameters [[stocks rs] dt]
  (let [start-params [0.1 0.1 0.1 0.1 0.1 0.1]
        proposal-cov  (ic/mult 0.001 (ic/identity-matrix 6)) 
        proposal-fn   (fn [old rand]
                        (dist/to-multivar-normal rand :mu (ic/matrix old)
                                                 :cov proposal-cov))
        proposal-dens (fn [old cand]
                        (Math/log
                         (dist/pdf-multivar-normal cand :mu (ic/matrix old)
                                                   :cov proposal-cov)))
        true-dens     (fn [x] (ll/black-scholes-vasicek-loglik x stocks rs dt))]
    (mcmc start-params proposal-fn proposal-dens true-dens
          (rand-stream) :logs true)))

(defn curved-propose [[s x] [rnd1 rnd2] & {:keys [base sd-c sd-x]
                                            :or {base 1.0 sd-c 1.0 sd-x 1.0}}]
   (let [c     (* s (Math/pow base x))
         new-c (dist/to-normal rnd1 :mean c :sd sd-c)
         new-x (dist/to-normal rnd2 :mean x :sd sd-x) 
         new-s (/ new-c (Math/pow base new-x))]
     [new-s new-x]))

(defn curved-dens
   [[old-s old-x] [new-s new-x] & {:keys [base sd-c sd-x] 
                                   :or {base 1.0 sd-c 1.0 sd-x 1.0}}]
   (let [old-c (* old-s (Math/pow base old-x))
         new-c (* new-s (Math/pow base new-x))]
     (* (is/pdf-normal new-c :mean old-c :sd sd-c)
        (is/pdf-normal new-x :mean old-x :sd sd-x))))

(defn cev-proposal-fn
  [avg-S [sd-mu sd-cS sd-al] [old-mu old-sS old-al] [rnd-mu rnd-sS rnd-al]]
  (concat
   [(dist/to-normal rnd-mu :mean old-mu :sd sd-mu)]
   (curved-propose [old-sS old-al] [rnd-sS rnd-al]
                   :base avg-S :sd-c sd-cS :sd-x sd-al)))

(defn cev-proposal-logdens
  [avg-S [sd-mu sd-cS sd-al] [old-mu old-sS old-al] [new-mu new-sS new-al]]
  (+ (Math/log (is/pdf-normal new-mu :mean old-mu :sd sd-mu))
     (Math/log (curved-dens [old-sS old-al] [new-sS new-al]
                                :base avg-S :sd-c sd-cS :sd-x sd-al))))

(defn ckls-proposal-fn
  [avg-r [sd-ka sd-be sd-cr sd-xi]
   [old-ka old-be old-sr old-xi] [rnd-ka rnd-be rnd-sr rnd-xi]]
  (concat
   [(dist/to-normal rnd-ka :mean old-ka :sd sd-ka)
    (dist/to-normal rnd-be :mean old-be :sd sd-be)]
   (curved-propose [old-sr old-xi] [rnd-sr rnd-xi]
                   :base avg-r :sd-c sd-cr :sd-x sd-xi)))

(defn ckls-proposal-logdens
  [avg-r [sd-ka sd-be sd-cr sd-xi]
   [old-ka old-be old-sr old-xi] [new-ka new-be new-sr new-xi]]
  (+ (Math/log (is/pdf-normal new-ka :mean old-ka :sd sd-ka))
     (Math/log (is/pdf-normal new-be :mean old-be :sd sd-be))
     (Math/log (curved-dens [old-sr old-xi] [new-sr new-xi]
                            :base avg-r :sd-c sd-cr :sd-x sd-xi))))

(defn cev-ckls-proposal-fn [avgs sds old rnd]
  (concat
   (cev-proposal-fn  (avgs 0) (! sds 0 3) (! old 0 3) (! rnd 0 3))  
   (ckls-proposal-fn (avgs 1) (! sds 3 7) (! old 3 7) (! rnd 3 7))
   [(dist/to-normal (nth rnd 7) :mean (nth old 7)  :sd (nth sds 7))]))

(defn cev-ckls-proposal-logdens [avgs sds old cnd]
  (+ (cev-proposal-logdens  (avgs 0) (! sds 0 3) (! old 0 3) (! cnd 0 3))
     (ckls-proposal-logdens (avgs 1) (! sds 3 7) (! old 3 7) (! cnd 3 7))
     (Math/log (is/pdf-normal (nth cnd 7) :mean (nth old 7) :sd (nth sds 7)))))

(defn sample-cev-parameters [stocks dt]
   (let [start-params  [0.05 0.05 1.0]
         sd-mu         0.01
         [sd-cS sd-al] [0.002 0.05]
         avg-S         (average stocks) 
         proposal-fn   (partial cev-proposal-fn avg-S [sd-mu sd-cS sd-al])
         proposal-dens (partial cev-proposal-logdens avg-S [sd-mu sd-cS sd-al])
         true-dens     (fn [x] (ll/cev-loglik x stocks dt))]
     (mcmc start-params proposal-fn proposal-dens true-dens (rand-stream) :logs true)))

(defn sample-ckls-parameters [rs dt]
  (let [start-params  [0.1 0.05 0.05 1.0]
        [sd-ka sd-be] [0.01 0.002]
        [sd-cr sd-xi] [0.001 0.02]
        avg-r         (average rs) 
        proposal-fn   (partial ckls-proposal-fn avg-r [sd-ka sd-be sd-cr sd-xi])
        proposal-dens (partial ckls-proposal-logdens avg-r [sd-ka sd-be sd-cr sd-xi])
        true-dens     (fn [x] (ll/ckls-loglik x rs dt))]
    (mcmc start-params proposal-fn proposal-dens true-dens (rand-stream) :logs true)))

(defn sample-cev-ckls-parameters [[stocks rs] dt]
  (let [start-params  [0.05 0.05 1.0 0.1 0.05 0.05 1.0 0.0]
        [sd-mu sd-ka sd-be sd-ro] [0.01 0.005 0.002 0.01]
        [sd-cS sd-cr sd-al sd-xi] [0.002 0.0005 0.05 0.01]
        [avg-S avg-r] [(average stocks) (average rs)]
        proposal-fn   (partial cev-ckls-proposal-fn [avg-S avg-r]
                              [sd-mu sd-cS sd-al sd-ka sd-be sd-cr sd-xi sd-ro])
        proposal-dens (partial cev-ckls-proposal-logdens [avg-S avg-r]
                              [sd-mu sd-cS sd-al sd-ka sd-be sd-cr sd-xi sd-ro])
        true-dens     (fn [x] (ll/cev-ckls-loglik x stocks rs dt))]
    (mcmc start-params proposal-fn proposal-dens true-dens
          (rand-stream) :logs true)))
