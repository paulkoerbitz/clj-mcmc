(ns clj-mcmc.helpers
  "Provides helper functions and macros")

(defmacro dbg [x]
  `(let [x# ~x] (println '~x "=" x#) x#))

(defn between [x low high]
  (min (max x low) high))
