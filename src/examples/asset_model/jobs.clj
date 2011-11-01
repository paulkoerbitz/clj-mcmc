((comment {:model   :bsv
   :dataset :1976-2006-monthly
   :outfile "results/2011-10-30/DaxMoneyCall_1976_2006_monthly_bs_vasicek_mcmc.csv"
   :n       10000
   :n-burn  5000
   :n-thin  10 }
  {:model   :bsv
   :dataset :1996-2006-monthly
   :outfile "results/2011-10-30/DaxMoneyCall_1996_2006_monthly_bs_vasicek_mcmc.csv"
   :n       10000
   :n-burn  5000
   :n-thin  10 }
  {:model   :bsv
   :dataset :1996-2006-weekly
   :outfile "results/2011-10-30/DaxMoneCall_1996_2006_weekly_bs_vasicek_mcmc.csv"
   :n       10000
   :n-burn  5000
   :n-thin  10 }
  {:model   :cev-ckls
   :dataset :1996-2006-monthly
   :outfile "results/2011-10-30/DaxMoneyCall_1996_2006_monthly_cev_ckls_mcmc.csv"
   :n       1000000
   :n-burn  5000
   :n-thin  10 }
  {:model   :cev-ckls
   :dataset :1976-2006-monthly
   :outfile "results/2011-10-30/DaxMoneyCall_1976_2006_monthly_cev_ckls_mcmc.csv"
   :n       1000000
   :n-burn  5000
   :n-thin  10 })
 {:model   :cev-ckls
  :dataset :1996-2006-weekly
  :outfile "results/2011-10-30/DaxMoneyCall_1996_2006_weekly_cev_ckls_mcmc.csv"
  :n       1000000
  :n-burn  5000
  :n-thin  10 }
 )
