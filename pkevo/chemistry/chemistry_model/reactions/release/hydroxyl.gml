rule [
  ruleID "Release from the synthase by nucleophilic attack of a hydroxyl group"
  left [
    edge [ source 2 target 3 label "-" ]
    edge [ source 6 target 7 label "-" ]
  ]
  context [
    node [ id 1 label "Ag" ]
    node [ id 2 label "S" ]
    node [ id 3 label "C" ]
    node [ id 4 label "O" ]
    edge [ source 1 target 2 label "-" ]
    edge [ source 3 target 4 label "=" ]
    node [ id 6 label "O" ]
    node [ id 7 label "H" ]
  ]
  right [
    edge [ source 2 target 7 label "-" ]
    edge [ source 6 target 3 label "-" ]
  ]
]
