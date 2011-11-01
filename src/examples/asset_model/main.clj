(ns examples.asset-model.main
  (:require [clojure-csv.core :as csv]
            [incanter.stats :as is]
            [incanter.core :as ic]
            [incanter.charts :as chrt]
            [clojure.java.io :as io])
  (:use [clj-mcmc core random distributions helpers]
        [examples.asset-model.mcmc-samplers]))

(set! *warn-on-reflection* true)

(def *dataset-file*
  "src/examples/asset_model/datasets.clj")

(def *job-file*
  "src/examples/asset_model/jobs.clj")

(defn parse-path-file [filename]
  (let [csv (csv/parse-csv (slurp filename))]
    [(float-array (for [[s r] csv :while r] (read-string s)))
     (float-array (for [[_ r] csv :while r] (read-string r)))]))

(defn show-scatter [coll i1 i2]
  (ic/view (chrt/scatter-plot (map #(nth % i1) coll) (map #(nth % i2) coll))))

(defn show-chart [chart-f coll i1]
  (ic/view (chart-f (map #(nth % i1) coll))))

(def show-hist (partial show-chart #(chrt/histogram % :nbins 100)))
(def show-trace (partial show-chart #(chrt/trace-plot %)))

(defn estimate-parameters
  ""
  [{:keys [n n-burn n-thin dataset model outfile]}]
  (let [dataset   (dataset (read-string (slurp *dataset-file*)))
        path      (parse-path-file (:location dataset))
        sample-fn (cond
                   (= model :cev)      sample-cev-parameters
                   (= model :ckls)     sample-ckls-parameters
                   (= model :bsv)      sample-bsv-parameters
                   (= model :cev-ckls) sample-cev-ckls-parameters
                   :else (throw (IllegalArgumentException.
                                 (str "Illegal model, valid models are "
                                      ":cev :ckls :bsv :cev-ckls."))))
        result     (take n (thin n-thin (drop n-burn (sample-fn path (:dt dataset)))))]
    (with-open [wrtr (io/writer outfile)]
      (.write wrtr "mu,ss,al,ka,be,sr,xi,ro")
      (doseq [line result]
        (.write wrtr
         (apply str (conj (interpose "," (map #(format "%9.7f" %) line)) "\n")))))
    result))

(defn -main []
  (let [jobs (read-string (slurp *job-file*))]
    (doseq [job jobs]
      (println "Doing job " job)
      (when (map? job) 
        (estimate-parameters job)))))
