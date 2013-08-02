(TeX-add-style-hook "introduction"
 (lambda ()
    (LaTeX-add-bibitems
     "Dutta2009"
     "Redestig2011")
    (LaTeX-add-labels
     "fig:box"
     "fig:bar"
     "fig:line")
    (TeX-run-style-hooks
     "a4wide"
     "hyperref"
     "latex2e"
     "art10"
     "article"
     "a4paper")))

