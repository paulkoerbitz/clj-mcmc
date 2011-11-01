(ns clj-mcmc.helpers
  "Provides helper functions and macros")

(defmacro dbg [x]
  `(let [x# ~x] (println '~x "=" x#) x#))

(defn between [x low high]
  (min (max x low) high))

(defn map-or-invoke [f coll]
  (if (coll? coll)
    (map f coll)
    (f coll)))

(defn average [xs]
  (/ (reduce + xs) (count xs)))

(defmacro prefix [pref items]
  (map-or-invoke #(symbol (str pref %)) items))

(defn nthcol 
  "Returns the nth entry of every item in coll.
   Example:
     (nthcol [[1 2 3] [4 5 6]] 1)
     ;=> (2 5)
  "
  [coll n]
  (for [x coll] (nth x n)))

(defmacro !
   "Returns a subvector [(v s) (v (+ s 1)) .. (v (- e 1))]
   Example:
     (! v 0 3)
     ;=> [(v 0) (v 1) (v 2)]
  "
   [v s e]
   (vec (for [i (range s e)] `(nth ~v ~i))))
