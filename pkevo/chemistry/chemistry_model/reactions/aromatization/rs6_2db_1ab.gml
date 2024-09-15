rule [
  ruleID "Aromatization  of a 6-membered ring with 2 double bonds and 1 aromatic bond"

  left [
    edge [source 1 target 2 label "-"]
    edge [source 2 target 3 label "="]
    edge [source 3 target 4 label "-"]
    edge [source 5 target 6 label "-"]
    edge [source 6 target 1 label "="]
  ]
  
  context [
    node [id 1 label "C"]
    node [id 2 label "C"]
    node [id 3 label "C"]
    node [id 6 label "C"]
    node [id 4 label "C"]
    node [id 5 label "C"]
    edge [source 4 target 5 label ":"]
  ]

  right [
    edge [source 1 target 2 label ":"]
    edge [source 2 target 3 label ":"]
    edge [source 3 target 4 label ":"]
    edge [source 5 target 6 label ":"]
    edge [source 6 target 1 label ":"]
  ]
]
