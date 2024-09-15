rule [
  ruleID "Aliphatic double bond oxidation"
  left [
    node [id 7 label "H+"]
    node [id 8 label "H-"]
    edge [source 3 target 4 label "-"]
    edge [source 5 target 6 label "="]
  ]
  context [
    node [id 0 label "C"]
    node [id 1 label "C"]
    node [id 2 label "H"]
    node [id 3 label "C"]
    node [id 4 label "H"]
    node [id 5 label "O"]
    node [id 6 label "O"]
    edge [source 0 target 3 label "-"]
    edge [source 1 target 3 label "-"]
    edge [source 2 target 3 label "-"]
  ]
  right [
    node [id 7 label "H"]
    node [id 8 label "H"]
    edge [source 3 target 5 label "-"]
    edge [source 5 target 4 label "-"]
    edge [source 6 target 7 label "-"]
    edge [source 6 target 8 label "-"]
  ]
]
