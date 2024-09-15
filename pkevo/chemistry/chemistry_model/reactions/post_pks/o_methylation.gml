rule [
  ruleID "O-Methylation"
  left [
    node [id 3 label "H"]
    node [id 4 label "C+"]
    edge [source 2 target 3 label "-"]
  ]
  context [
    node [id 1 label "C"]
    node [id 2 label "O"]
    node [id 5 label "H"]
    node [id 6 label "H"]
    node [id 7 label "H"]
    edge [source 1 target 2 label "-"]
    edge [source 4 target 5 label "-"]
    edge [source 4 target 6 label "-"]
    edge [source 4 target 7 label "-"]
  ]
  right [
    node [id 3 label "H+"]
    node [id 4 label "C"]
    edge [source 2 target 4 label "-"]
  ]
]
