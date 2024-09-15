rule [
  ruleID "Aromatization  of a 6-membered ring with 1 carbonyl group, 1 double bond and 1 aromatic bond - orientation 3"

  left [
    edge [source 1 target 2 label "-"]
    edge [source 3 target 4 label "-"]
    edge [source 4 target 5 label "-"]
    edge [source 5 target 6 label "="]
    edge [source 6 target 1 label "-"]

    edge [source 1 target 7 label "="]
    edge [source 4 target 8 label "-"]
  ]
  
  context [
    node [id 1 label "C"]
    node [id 4 label "C"]
    node [id 5 label "C"]
    node [id 6 label "C"]
    node [id 2 label "C"]
    node [id 3 label "C"]
    edge [source 2 target 3 label ":"]

    node [id 7 label "O"]
    node [id 8 label "H"]
  ]

  right [
    edge [source 1 target 2 label ":"]
    edge [source 3 target 4 label ":"]
    edge [source 4 target 5 label ":"]
    edge [source 5 target 6 label ":"]
    edge [source 6 target 1 label ":"]

    edge [source 1 target 7 label "-"]
    edge [source 7 target 8 label "-"]
  ]
]
