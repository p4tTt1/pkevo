rule [
  ruleID "Intramolecular Claisen thioester condensation via near carbonyl group, ringsize 6"

  left [
    edge [source 1 target 11 label "-"]
    edge [source 6 target 10 label "-"]
  ]
  
  context [
    node [id  1 label "C"]
    node [id  2 label "*"]
    node [id  3 label "*"]
    node [id  4 label "*"]
    node [id  5 label "C"]
    node [id  6 label "C"]
    node [id  7 label "C"]
    edge [source 1 target 2 label "*"]
    edge [source 2 target 3 label "*"]
    edge [source 3 target 4 label "*"]
    edge [source 4 target 5 label "*"]
    edge [source 5 target 6 label "-"]
    edge [source 6 target 7 label "-"]

    node [id  8 label "O"]
    node [id  9 label "O"]
    edge [source 1 target 8 label "="]
    edge [source 7 target 9 label "="]

    node [id 10 label "H"]
    node [id 11 label "S"]
  ]

  right [
    edge [source 1 target 6 label "-"]
    edge [source 11 target 10 label "-"]
  ]
  
  constrainShortestPath [
    source 1 target 6
    op "=" length 5
  ]
]
