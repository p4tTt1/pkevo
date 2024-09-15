rule [
  ruleID "Intramolecular aldol addition - ringsize 6 - distant carbonyl group"

  left [
    edge [source 7 target 9 label "="]
    edge [source 2 target 10 label "-"]
  ]

  context [
    node [id  1 label "C"]
    node [id  2 label "C"]
    node [id  3 label "*"]
    node [id  4 label "*"]
    node [id  5 label "*"]
    node [id  6 label "*"]
    node [id  7 label "C"]
    node [id  8 label "O"]
    node [id  9 label "O"]
    node [id 10 label "H"]
    edge [source 1 target 2 label "-"]
    edge [source 2 target 3 label "*"]
    edge [source 3 target 4 label "*"]
    edge [source 4 target 5 label "*"]
    edge [source 5 target 6 label "*"]
    edge [source 6 target 7 label "*"]

    edge [source 1 target 8 label "="]
  ]

  right [
    edge [source  2 target  7 label "-"]
    edge [source  7 target  9 label "-"]
    edge [source  9 target  10 label "-"]
  ]

  constrainShortestPath [
    source 2 target 7
    op "=" length 5
  ]
]
