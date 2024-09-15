rule [
  ruleID "Epoxidation"
  left [
    node [id 5 label "H+"]
    node [id 6 label "H-"]
    edge [source 1 target 2 label "="]
    edge [source 3 target 4 label "="]
  ]
  context [
    node [id 1 label "C"]
    node [id 2 label "C"]
    node [id 3 label "O"]
    node [id 4 label "O"]
  ]
  right [
    node [id 5 label "H"]
    node [id 6 label "H"]
    edge [source 1 target 2 label "-"]
    edge [source 1 target 3 label "-"]
    edge [source 2 target 3 label "-"]
    edge [source 4 target 5 label "-"]
    edge [source 4 target 6 label "-"]
  ]
]
