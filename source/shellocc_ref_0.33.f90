!! ------------------------------------------------------------------------
!! neutral atomic shell occupations at wB97X-V/vDZP
!! average over fit set molecules
!! gtb coord -apo -nogtb
!! MLmix = 1/3 
!! ------------------------------------------------------------------------

subroutine shellocc_ref(at,socc) 
      use bascom
      implicit none
      integer at
      real*8 socc(*)

      socc(1:bas_nsh(at))=0

      select case (at)
      case (1)
      socc( 1)=  0.69264222808434
      socc( 2)=  0.14159955297977
      socc( 3)=  0.16575821893589
      case (2)
      socc( 1)=  1.30320944384752
      socc( 2)=  0.67445566751767
      socc( 3)=  0.02233488863481
      case (3)
      socc( 1)=  1.95899377872121
      socc( 2)=  0.25267315485192
      socc( 3)=  0.11620627212046
      socc( 4)=  0.41258386001878
      socc( 5)=  0.25954293428763
      case (4)
      socc( 1)=  1.97111378215842
      socc( 2)=  0.45810747906670
      socc( 3)=  0.15300493595427
      socc( 4)=  0.94930693774026
      socc( 5)=  0.46846686508036
      case (5)
      socc( 1)=  0.48354388386042
      socc( 2)=  0.15681644999964
      socc( 3)=  1.47322175838681
      socc( 4)=  0.38482788015938
      socc( 5)=  0.50159002759375
      case (6)
      socc( 1)=  0.55666397324469
      socc( 2)=  0.25120261811235
      socc( 3)=  2.10755997506139
      socc( 4)=  0.75746522555791
      socc( 5)=  0.32710820802367
      case (7)
      socc( 1)=  0.79624139865659
      socc( 2)=  0.31749496440959
      socc( 3)=  2.73635251741813
      socc( 4)=  0.98096348058478
      socc( 5)=  0.16894763893091
      case (8)
      socc( 1)=  1.08460674785680
      socc( 2)=  0.38134064737674
      socc( 3)=  3.32762843994261
      socc( 4)=  1.12937855100283
      socc( 5)=  0.07704561382102
      case (9)
      socc( 1)=  1.26298582965783
      socc( 2)=  0.45839005699390
      socc( 3)=  3.83052034186561
      socc( 4)=  1.40814828911248
      socc( 5)=  0.03995548237017
      case (10)
      socc( 1)=  1.47024479408269
      socc( 2)=  0.52146835610379
      socc( 3)=  3.70575591161481
      socc( 4)=  2.28623261801831
      socc( 5)=  0.01629832018040
      case (11)
      socc( 1)=  1.44272104623990
      socc( 2)=  0.56147696803778
      socc( 3)=  0.30286795060414
      socc( 4)=  5.15471618358905
      socc( 5)=  1.23531381129067
      socc( 6)=  0.30290404023845
      case (12)
      socc( 1)=  1.41118194136168
      socc( 2)=  0.60257450561694
      socc( 3)=  0.58202476396721
      socc( 4)=  5.90203159422521
      socc( 5)=  0.88913520405121
      socc( 6)=  0.61305199077774
      case (13)
      socc( 1)=  0.55409123445529
      socc( 2)=  0.18700939319552
      socc( 3)=  1.14200108985144
      socc( 4)=  0.36773005186021
      socc( 5)=  0.74916823063754
      case (14)
      socc( 1)=  0.78489160465251
      socc( 2)=  0.22437622129695
      socc( 3)=  1.89887560555902
      socc( 4)=  0.54518793303536
      socc( 5)=  0.54666863545615
      case (15)
      socc( 1)=  0.99441180454863
      socc( 2)=  0.38348294255053
      socc( 3)=  2.42166885021688
      socc( 4)=  0.62108522469704
      socc( 5)=  0.57935117798692
      case (16)
      socc( 1)=  1.16875462402413
      socc( 2)=  0.42250722184002
      socc( 3)=  3.18695750632506
      socc( 4)=  0.93092404067175
      socc( 5)=  0.29085660713904
      case (17)
      socc( 1)=  1.26904700169420
      socc( 2)=  0.50124825150823
      socc( 3)=  3.96697497126814
      socc( 4)=  1.15243762904417
      socc( 5)=  0.11029214648526
      case (18)
      socc( 1)=  1.36124497167220
      socc( 2)=  0.63500946581722
      socc( 3)=  4.81587227054938
      socc( 4)=  1.06582014894422
      socc( 5)=  0.12205314301698
      case (19) 
      case (20) 
      case (21) 
      case (22)
      case (31) 
      case (32) 
      case (33) 
      case (34) 
      case (35) 
      socc( 1)=  1.40434899781744
      socc( 2)=  0.42022034962529
      socc( 3)=  3.94291335821524
      socc( 4)=  1.12941461820133
      socc( 5)=  0.10310267614070
      case (36) 
      case (37) 
      case (38) 
      case (46) 
      case (49) 
      case (50) 
      case (51)
      case (52) 
      case (53) 
      socc( 1)=  1.24287604544147
      socc( 2)=  0.61429622795139
      socc( 3)=  3.93864154929533
      socc( 4)=  1.08527271224683
      socc( 5)=  0.11891346506498
      case (54)
      case (55) 
      case (56) 
      case (81)
      case (82)
      case (83)
      case (84)
      case (85)
      case (86) 
      end select

      end

