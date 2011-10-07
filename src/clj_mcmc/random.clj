(ns clj-mcmc.random
  "clj-mcmc.random warps various random samplers as lazy-seqs
   to simplify functional usage in sampling algorithms and
   separation of concerns."
  )

(defn rand-stream []
  (cons (rand) (lazy-seq (rand-stream))))
