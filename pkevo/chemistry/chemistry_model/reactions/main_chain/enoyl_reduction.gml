rule [
  ruleID "Reduction of enoyl group"
  left [
    node [ id 7 label "H-" ]
    node [ id 8 label "H+" ]
    edge [ source 4 target 5 label "=" ]
  ]
  context [
    node [ id 0 label "Ag" ]
    node [ id 1 label "S" ]
    node [ id 2 label "C" ]
    node [ id 3 label "O" ]
    node [ id 4 label "C" ]
    node [ id 5 label "C" ]
    node [ id 6 label "C" ]
    edge [ source 0 target 1 label "-" ]
    edge [ source 1 target 2 label "-" ]
    edge [ source 2 target 3 label "=" ]
    edge [ source 2 target 4 label "-" ]
    edge [ source 5 target 6 label "-" ]
  ]
  right [
    node [ id 7 label "H" ]
    node [ id 8 label "H" ]
    edge [ source 4 target 5 label "-" ]
    edge [ source 4 target 8 label "-" ]
    edge [ source 5 target 7 label "-" ]
  ]
]
