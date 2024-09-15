rule [
  ruleID "Aromatization  of a 6-membered ring with 2 carbonyl groups and 1 hydroxyl group"

  left [
    edge [source 1 target 2 label "-"]
    edge [source 2 target 3 label "-"]
    edge [source 3 target 4 label "-"]
    edge [source 4 target 5 label "-"]
    edge [source 5 target 6 label "-"]
    edge [source 6 target 1 label "-"]

    edge [source 1 target 7 label "="]
    edge [source 2 target 10 label "-"]
    edge [source 3 target 8 label "="]
    edge [source 4 target 11 label "-"]
    edge [source 5 target 9 label "-"]
    edge [source 6 target 12 label "-"]
  ]
  
  context [
    node [id 1 label "C"]
    node [id 2 label "C"]
    node [id 3 label "C"]
    node [id 4 label "C"]
    node [id 5 label "C"]
    node [id 6 label "C"]
    node [id 7 label "O"]
    node [id 8 label "O"]
    node [id 9 label "O"]
    node [id 10 label "H"]
    node [id 11 label "H"]
    node [id 12 label "H"]
    node [id 13 label "H"]
    edge [source 9 target 13 label "-"]
  ]

  right [
    edge [source 1 target 2 label ":"]
    edge [source 2 target 3 label ":"]
    edge [source 3 target 4 label ":"]
    edge [source 4 target 5 label ":"]
    edge [source 5 target 6 label ":"]
    edge [source 6 target 1 label ":"]

    edge [source 1 target 7 label "-"]
    edge [source 7 target 10 label "-"]
    edge [source 3 target 8 label "-"]
    edge [source 8 target 11 label "-"]
    edge [source 9 target 12 label "-"]
  ]
]
