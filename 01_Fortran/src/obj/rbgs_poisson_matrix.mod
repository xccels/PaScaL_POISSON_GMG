V34 :0x24 rbgs_poisson_matrix
23 rbgs_poisson_matrix.f90 S624 0
05/04/2026  18:30:16
use geometry private
use matrix private
enduse
D 82 26 674 1448 673 7
D 145 22 7
D 147 22 7
D 149 22 7
D 151 22 7
D 153 22 7
D 155 22 7
D 157 22 7
D 159 22 7
D 161 22 7
D 166 26 741 2288 740 7
D 250 22 7
D 252 22 7
D 254 22 7
D 256 22 7
D 258 22 7
D 260 22 7
D 262 22 7
D 264 22 7
D 266 22 7
D 268 22 7
D 270 22 7
D 272 22 7
D 322 26 1070 304 1069 7
D 334 22 7
D 491 23 10 3 445 457 1 1 0 0 1
 10 446 11 447 446 448
 10 449 450 451 449 452
 10 453 454 455 453 456
D 494 23 10 3 458 470 1 1 0 0 1
 10 459 11 460 459 461
 10 462 463 464 462 465
 10 466 467 468 466 469
D 497 23 18 1 10 13 0 0 0 0 0
 10 12 11 10 12 13
D 500 23 10 3 471 483 1 1 0 0 1
 10 472 11 473 472 474
 10 475 476 477 475 478
 10 479 480 481 479 482
D 503 23 10 3 484 496 1 1 0 0 1
 10 485 11 486 485 487
 10 488 489 490 488 491
 10 492 493 494 492 495
D 506 23 18 1 10 13 0 0 0 0 0
 10 12 11 10 12 13
S 624 24 0 0 0 9 1 0 5013 10015 0 A 0 0 0 0 B 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 rbgs_poisson_matrix
S 627 3 0 0 0 7 1 1 0 0 0 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 2 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 7
S 628 3 0 0 0 7 1 1 0 0 0 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 3 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 7
S 661 3 0 0 0 7 1 1 0 0 0 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 17 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 7
S 667 3 0 0 0 7 1 1 0 0 0 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 29 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 7
R 673 25 1 geometry domain
R 674 5 2 geometry nx domain
R 675 5 3 geometry ny domain
R 676 5 4 geometry nz domain
R 677 5 5 geometry lx domain
R 678 5 6 geometry ly domain
R 679 5 7 geometry lz domain
R 680 5 8 geometry ox domain
R 681 5 9 geometry oy domain
R 682 5 10 geometry oz domain
R 684 5 12 geometry dxm domain
R 685 5 13 geometry dxm$sd domain
R 686 5 14 geometry dxm$p domain
R 687 5 15 geometry dxm$o domain
R 689 5 17 geometry dym domain
R 691 5 19 geometry dym$sd domain
R 692 5 20 geometry dym$p domain
R 693 5 21 geometry dym$o domain
R 695 5 23 geometry dzm domain
R 697 5 25 geometry dzm$sd domain
R 698 5 26 geometry dzm$p domain
R 699 5 27 geometry dzm$o domain
R 702 5 30 geometry dxg domain
R 703 5 31 geometry dxg$sd domain
R 704 5 32 geometry dxg$p domain
R 705 5 33 geometry dxg$o domain
R 707 5 35 geometry dyg domain
R 709 5 37 geometry dyg$sd domain
R 710 5 38 geometry dyg$p domain
R 711 5 39 geometry dyg$o domain
R 713 5 41 geometry dzg domain
R 715 5 43 geometry dzg$sd domain
R 716 5 44 geometry dzg$p domain
R 717 5 45 geometry dzg$o domain
R 720 5 48 geometry xg domain
R 721 5 49 geometry xg$sd domain
R 722 5 50 geometry xg$p domain
R 723 5 51 geometry xg$o domain
R 725 5 53 geometry yg domain
R 727 5 55 geometry yg$sd domain
R 728 5 56 geometry yg$p domain
R 729 5 57 geometry yg$o domain
R 731 5 59 geometry zg domain
R 733 5 61 geometry zg$sd domain
R 734 5 62 geometry zg$p domain
R 735 5 63 geometry zg$o domain
R 737 5 65 geometry is_periodic domain
R 740 25 68 geometry subdomain
R 741 5 69 geometry nx subdomain
R 742 5 70 geometry ny subdomain
R 743 5 71 geometry nz subdomain
R 744 5 72 geometry lx subdomain
R 745 5 73 geometry ly subdomain
R 746 5 74 geometry lz subdomain
R 747 5 75 geometry ox subdomain
R 748 5 76 geometry oy subdomain
R 749 5 77 geometry oz subdomain
R 751 5 79 geometry dxm subdomain
R 752 5 80 geometry dxm$sd subdomain
R 753 5 81 geometry dxm$p subdomain
R 754 5 82 geometry dxm$o subdomain
R 757 5 85 geometry dym subdomain
R 758 5 86 geometry dym$sd subdomain
R 759 5 87 geometry dym$p subdomain
R 760 5 88 geometry dym$o subdomain
R 763 5 91 geometry dzm subdomain
R 764 5 92 geometry dzm$sd subdomain
R 765 5 93 geometry dzm$p subdomain
R 766 5 94 geometry dzm$o subdomain
R 769 5 97 geometry dxg subdomain
R 770 5 98 geometry dxg$sd subdomain
R 771 5 99 geometry dxg$p subdomain
R 772 5 100 geometry dxg$o subdomain
R 775 5 103 geometry dyg subdomain
R 776 5 104 geometry dyg$sd subdomain
R 777 5 105 geometry dyg$p subdomain
R 778 5 106 geometry dyg$o subdomain
R 781 5 109 geometry dzg subdomain
R 782 5 110 geometry dzg$sd subdomain
R 783 5 111 geometry dzg$p subdomain
R 784 5 112 geometry dzg$o subdomain
R 787 5 115 geometry xg subdomain
R 788 5 116 geometry xg$sd subdomain
R 789 5 117 geometry xg$p subdomain
R 790 5 118 geometry xg$o subdomain
R 793 5 121 geometry yg subdomain
R 794 5 122 geometry yg$sd subdomain
R 795 5 123 geometry yg$p subdomain
R 796 5 124 geometry yg$o subdomain
R 799 5 127 geometry zg subdomain
R 800 5 128 geometry zg$sd subdomain
R 801 5 129 geometry zg$p subdomain
R 802 5 130 geometry zg$o subdomain
R 807 5 135 geometry x subdomain
R 808 5 136 geometry x$sd subdomain
R 809 5 137 geometry x$p subdomain
R 810 5 138 geometry x$o subdomain
R 812 5 140 geometry b subdomain
R 816 5 144 geometry b$sd subdomain
R 817 5 145 geometry b$p subdomain
R 818 5 146 geometry b$o subdomain
R 820 5 148 geometry r subdomain
R 824 5 152 geometry r$sd subdomain
R 825 5 153 geometry r$p subdomain
R 826 5 154 geometry r$o subdomain
R 828 5 156 geometry is_periodic subdomain
R 829 5 157 geometry is_aggregated subdomain
R 830 5 158 geometry ista subdomain
R 831 5 159 geometry iend subdomain
R 832 5 160 geometry jsta subdomain
R 833 5 161 geometry jend subdomain
R 834 5 162 geometry ksta subdomain
R 835 5 163 geometry kend subdomain
R 836 5 164 geometry ddt_yz_plane_x0 subdomain
R 837 5 165 geometry ddt_yz_plane_x1 subdomain
R 838 5 166 geometry ddt_yz_plane_xn subdomain
R 839 5 167 geometry ddt_yz_plane_xn1 subdomain
R 840 5 168 geometry ddt_xz_plane_y0 subdomain
R 841 5 169 geometry ddt_xz_plane_y1 subdomain
R 842 5 170 geometry ddt_xz_plane_yn subdomain
R 843 5 171 geometry ddt_xz_plane_yn1 subdomain
R 844 5 172 geometry ddt_xy_plane_z0 subdomain
R 845 5 173 geometry ddt_xy_plane_z1 subdomain
R 846 5 174 geometry ddt_xy_plane_zn subdomain
R 847 5 175 geometry ddt_xy_plane_zn1 subdomain
R 848 5 176 geometry ddt_inner_domain subdomain
R 849 5 177 geometry is_x0_boundary subdomain
R 850 5 178 geometry is_x1_boundary subdomain
R 851 5 179 geometry is_y0_boundary subdomain
R 852 5 180 geometry is_y1_boundary subdomain
R 853 5 181 geometry is_z0_boundary subdomain
R 854 5 182 geometry is_z1_boundary subdomain
S 1066 3 0 0 0 7 1 1 0 0 0 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 35 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 7
R 1069 25 1 matrix matrix_poisson
R 1070 5 2 matrix dof matrix_poisson
R 1071 5 3 matrix coeff matrix_poisson
R 1076 5 8 matrix coeff$sd matrix_poisson
R 1077 5 9 matrix coeff$p matrix_poisson
R 1078 5 10 matrix coeff$o matrix_poisson
S 1159 27 0 0 0 9 1161 624 7814 0 8000000 A 0 0 0 0 B 0 9 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 624 0 0 0 0 rbgs_solver_poisson_matrix
S 1160 27 0 0 0 9 1197 624 7841 0 8000000 A 0 0 0 0 B 0 10 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 624 0 0 0 0 rbgs_iterator_poisson_matrix
S 1161 23 5 0 0 0 1170 624 7814 0 0 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 rbgs_solver_poisson_matrix
S 1162 7 3 3 0 491 1 1161 7870 20000014 10003000 A 0 0 0 0 B 0 14 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 sol
S 1163 1 3 1 0 322 1 1161 7622 14 3000 A 0 0 0 0 B 0 14 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 a_poisson
S 1164 7 3 1 0 494 1 1161 7874 20000014 10003000 A 0 0 0 0 B 0 14 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 rhs
S 1165 1 3 1 0 166 1 1161 7757 14 3000 A 0 0 0 0 B 0 14 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 dm
S 1166 1 3 1 0 6 1 1161 7878 14 3000 A 0 0 0 0 B 0 14 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 maxiteration
S 1167 1 3 1 0 10 1 1161 7891 14 3000 A 0 0 0 0 B 0 14 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 tolerance
S 1168 1 3 1 0 10 1 1161 7901 14 3000 A 0 0 0 0 B 0 14 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 omega
S 1169 7 3 1 0 497 1 1161 6098 800014 3000 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 is_aggregated
S 1170 14 5 0 0 0 1 1161 7814 20000000 400000 A 0 0 0 0 B 0 0 0 0 0 0 0 75 8 0 0 0 0 0 0 0 0 0 0 0 0 14 0 624 0 0 0 0 rbgs_solver_poisson_matrix rbgs_solver_poisson_matrix 
F 1170 8 1162 1163 1164 1165 1166 1167 1168 1169
S 1171 6 1 0 0 7 1 1161 6708 40800016 3000 A 0 0 0 0 B 0 22 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 z_b_0_1
S 1172 6 1 0 0 7 1 1161 6716 40800016 3000 A 0 0 0 0 B 0 22 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 z_b_1_1
S 1173 6 1 0 0 7 1 1161 6724 40800016 3000 A 0 0 0 0 B 0 22 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 z_b_2_1
S 1174 6 1 0 0 7 1 1161 6732 40800016 3000 A 0 0 0 0 B 0 22 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 z_b_3_1
S 1175 6 1 0 0 7 1 1161 6740 40800016 3000 A 0 0 0 0 B 0 22 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 z_b_4_1
S 1176 6 1 0 0 7 1 1161 6748 40800016 3000 A 0 0 0 0 B 0 22 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 z_b_5_1
S 1177 6 1 0 0 7 1 1161 6756 40800016 3000 A 0 0 0 0 B 0 22 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 z_b_6_1
S 1178 6 1 0 0 7 1 1161 6764 40800016 3000 A 0 0 0 0 B 0 22 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 z_b_7_1
S 1179 6 1 0 0 7 1 1161 6772 40800016 3000 A 0 0 0 0 B 0 22 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 z_b_8_1
S 1180 6 1 0 0 7 1 1161 6780 40800016 3000 A 0 0 0 0 B 0 22 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 z_b_9_1
S 1181 6 1 0 0 7 1 1161 7907 40800016 3000 A 0 0 0 0 B 0 22 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 z_e_457
S 1182 6 1 0 0 7 1 1161 7915 40800016 3000 A 0 0 0 0 B 0 22 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 z_e_460
S 1183 6 1 0 0 7 1 1161 7923 40800016 3000 A 0 0 0 0 B 0 22 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 z_e_463
S 1184 6 1 0 0 7 1 1161 7025 40800016 3000 A 0 0 0 0 B 0 24 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 z_b_10_1
S 1185 6 1 0 0 7 1 1161 7034 40800016 3000 A 0 0 0 0 B 0 24 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 z_b_11_1
S 1186 6 1 0 0 7 1 1161 7043 40800016 3000 A 0 0 0 0 B 0 24 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 z_b_12_1
S 1187 6 1 0 0 7 1 1161 7052 40800016 3000 A 0 0 0 0 B 0 24 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 z_b_13_1
S 1188 6 1 0 0 7 1 1161 7061 40800016 3000 A 0 0 0 0 B 0 24 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 z_b_14_1
S 1189 6 1 0 0 7 1 1161 7070 40800016 3000 A 0 0 0 0 B 0 24 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 z_b_15_1
S 1190 6 1 0 0 7 1 1161 7079 40800016 3000 A 0 0 0 0 B 0 24 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 z_b_16_1
S 1191 6 1 0 0 7 1 1161 7088 40800016 3000 A 0 0 0 0 B 0 24 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 z_b_17_1
S 1192 6 1 0 0 7 1 1161 7097 40800016 3000 A 0 0 0 0 B 0 24 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 z_b_18_1
S 1193 6 1 0 0 7 1 1161 7106 40800016 3000 A 0 0 0 0 B 0 24 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 z_b_19_1
S 1194 6 1 0 0 7 1 1161 7931 40800016 3000 A 0 0 0 0 B 0 24 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 z_e_476
S 1195 6 1 0 0 7 1 1161 7939 40800016 3000 A 0 0 0 0 B 0 24 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 z_e_479
S 1196 6 1 0 0 7 1 1161 7947 40800016 3000 A 0 0 0 0 B 0 24 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 z_e_482
S 1197 23 5 0 0 0 1205 624 7841 0 0 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 rbgs_iterator_poisson_matrix
S 1198 7 3 3 0 500 1 1197 7870 20000014 10003000 A 0 0 0 0 B 0 119 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 sol
S 1199 1 3 1 0 322 1 1197 7622 14 3000 A 0 0 0 0 B 0 119 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 a_poisson
S 1200 7 3 1 0 503 1 1197 7874 20000014 10003000 A 0 0 0 0 B 0 119 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 rhs
S 1201 1 3 1 0 166 1 1197 7757 14 3000 A 0 0 0 0 B 0 119 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 dm
S 1202 1 3 1 0 6 1 1197 7878 14 3000 A 0 0 0 0 B 0 119 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 maxiteration
S 1203 1 3 1 0 10 1 1197 7901 14 3000 A 0 0 0 0 B 0 119 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 omega
S 1204 7 3 1 0 506 1 1197 6098 800014 3000 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 is_aggregated
S 1205 14 5 0 0 0 1 1197 7841 20000000 400000 A 0 0 0 0 B 0 0 0 0 0 0 0 84 7 0 0 0 0 0 0 0 0 0 0 0 0 119 0 624 0 0 0 0 rbgs_iterator_poisson_matrix rbgs_iterator_poisson_matrix 
F 1205 7 1198 1199 1200 1201 1202 1203 1204
S 1206 6 1 0 0 7 1 1197 6708 40800016 3000 A 0 0 0 0 B 0 127 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 z_b_0_1
S 1207 6 1 0 0 7 1 1197 6716 40800016 3000 A 0 0 0 0 B 0 127 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 z_b_1_1
S 1208 6 1 0 0 7 1 1197 6724 40800016 3000 A 0 0 0 0 B 0 127 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 z_b_2_1
S 1209 6 1 0 0 7 1 1197 6732 40800016 3000 A 0 0 0 0 B 0 127 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 z_b_3_1
S 1210 6 1 0 0 7 1 1197 6740 40800016 3000 A 0 0 0 0 B 0 127 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 z_b_4_1
S 1211 6 1 0 0 7 1 1197 6748 40800016 3000 A 0 0 0 0 B 0 127 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 z_b_5_1
S 1212 6 1 0 0 7 1 1197 6756 40800016 3000 A 0 0 0 0 B 0 127 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 z_b_6_1
S 1213 6 1 0 0 7 1 1197 6764 40800016 3000 A 0 0 0 0 B 0 127 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 z_b_7_1
S 1214 6 1 0 0 7 1 1197 6772 40800016 3000 A 0 0 0 0 B 0 127 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 z_b_8_1
S 1215 6 1 0 0 7 1 1197 6780 40800016 3000 A 0 0 0 0 B 0 127 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 z_b_9_1
S 1216 6 1 0 0 7 1 1197 7955 40800016 3000 A 0 0 0 0 B 0 127 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 z_e_483
S 1217 6 1 0 0 7 1 1197 7963 40800016 3000 A 0 0 0 0 B 0 127 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 z_e_486
S 1218 6 1 0 0 7 1 1197 7971 40800016 3000 A 0 0 0 0 B 0 127 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 z_e_489
S 1219 6 1 0 0 7 1 1197 7025 40800016 3000 A 0 0 0 0 B 0 129 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 z_b_10_1
S 1220 6 1 0 0 7 1 1197 7034 40800016 3000 A 0 0 0 0 B 0 129 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 z_b_11_1
S 1221 6 1 0 0 7 1 1197 7043 40800016 3000 A 0 0 0 0 B 0 129 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 z_b_12_1
S 1222 6 1 0 0 7 1 1197 7052 40800016 3000 A 0 0 0 0 B 0 129 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 z_b_13_1
S 1223 6 1 0 0 7 1 1197 7061 40800016 3000 A 0 0 0 0 B 0 129 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 z_b_14_1
S 1224 6 1 0 0 7 1 1197 7070 40800016 3000 A 0 0 0 0 B 0 129 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 z_b_15_1
S 1225 6 1 0 0 7 1 1197 7079 40800016 3000 A 0 0 0 0 B 0 129 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 z_b_16_1
S 1226 6 1 0 0 7 1 1197 7088 40800016 3000 A 0 0 0 0 B 0 129 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 z_b_17_1
S 1227 6 1 0 0 7 1 1197 7097 40800016 3000 A 0 0 0 0 B 0 129 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 z_b_18_1
S 1228 6 1 0 0 7 1 1197 7106 40800016 3000 A 0 0 0 0 B 0 129 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 z_b_19_1
S 1229 6 1 0 0 7 1 1197 7979 40800016 3000 A 0 0 0 0 B 0 129 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 z_e_502
S 1230 6 1 0 0 7 1 1197 7987 40800016 3000 A 0 0 0 0 B 0 129 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 z_e_505
S 1231 6 1 0 0 7 1 1197 7995 40800016 3000 A 0 0 0 0 B 0 129 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 z_e_508
A 12 2 0 0 0 7 627 0 0 0 12 0 0 0 0 0 0 0 0 0 0 0
A 13 2 0 0 0 7 628 0 0 0 13 0 0 0 0 0 0 0 0 0 0 0
A 14 2 0 0 0 7 661 0 0 0 14 0 0 0 0 0 0 0 0 0 0 0
A 164 2 0 0 0 7 667 0 0 0 164 0 0 0 0 0 0 0 0 0 0 0
A 367 2 0 0 0 7 1066 0 0 0 367 0 0 0 0 0 0 0 0 0 0 0
A 445 1 0 0 0 7 1180 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
A 446 1 0 0 0 7 1172 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
A 447 1 0 0 0 7 1171 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
A 448 1 0 0 0 7 1181 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
A 449 1 0 0 0 7 1175 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
A 450 1 0 0 0 7 1173 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
A 451 1 0 0 0 7 1174 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
A 452 1 0 0 0 7 1182 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
A 453 1 0 0 0 7 1178 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
A 454 1 0 0 0 7 1176 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
A 455 1 0 0 0 7 1177 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
A 456 1 0 0 0 7 1183 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
A 457 1 0 0 0 7 1179 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
A 458 1 0 0 0 7 1193 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
A 459 1 0 0 0 7 1185 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
A 460 1 0 0 0 7 1184 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
A 461 1 0 0 0 7 1194 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
A 462 1 0 0 0 7 1188 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
A 463 1 0 0 0 7 1186 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
A 464 1 0 0 0 7 1187 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
A 465 1 0 0 0 7 1195 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
A 466 1 0 0 0 7 1191 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
A 467 1 0 0 0 7 1189 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
A 468 1 0 0 0 7 1190 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
A 469 1 0 0 0 7 1196 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
A 470 1 0 0 0 7 1192 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
A 471 1 0 0 0 7 1215 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
A 472 1 0 0 0 7 1207 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
A 473 1 0 0 0 7 1206 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
A 474 1 0 0 0 7 1216 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
A 475 1 0 0 0 7 1210 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
A 476 1 0 0 0 7 1208 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
A 477 1 0 0 0 7 1209 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
A 478 1 0 0 0 7 1217 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
A 479 1 0 0 0 7 1213 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
A 480 1 0 0 0 7 1211 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
A 481 1 0 0 0 7 1212 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
A 482 1 0 0 0 7 1218 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
A 483 1 0 0 0 7 1214 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
A 484 1 0 0 0 7 1228 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
A 485 1 0 0 0 7 1220 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
A 486 1 0 0 0 7 1219 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
A 487 1 0 0 0 7 1229 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
A 488 1 0 0 0 7 1223 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
A 489 1 0 0 0 7 1221 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
A 490 1 0 0 0 7 1222 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
A 491 1 0 0 0 7 1230 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
A 492 1 0 0 0 7 1226 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
A 493 1 0 0 0 7 1224 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
A 494 1 0 0 0 7 1225 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
A 495 1 0 0 0 7 1231 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
A 496 1 0 0 0 7 1227 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
Z
T 673 82 0 0 0 0
A 686 7 145 0 1 2 1
A 685 7 0 14 1 10 1
A 692 7 147 0 1 2 1
A 691 7 0 14 1 10 1
A 698 7 149 0 1 2 1
A 697 7 0 14 1 10 1
A 704 7 151 0 1 2 1
A 703 7 0 14 1 10 1
A 710 7 153 0 1 2 1
A 709 7 0 14 1 10 1
A 716 7 155 0 1 2 1
A 715 7 0 14 1 10 1
A 722 7 157 0 1 2 1
A 721 7 0 14 1 10 1
A 728 7 159 0 1 2 1
A 727 7 0 14 1 10 1
A 734 7 161 0 1 2 1
A 733 7 0 14 1 10 0
T 740 166 0 0 0 0
A 753 7 250 0 1 2 1
A 752 7 0 14 1 10 1
A 759 7 252 0 1 2 1
A 758 7 0 14 1 10 1
A 765 7 254 0 1 2 1
A 764 7 0 14 1 10 1
A 771 7 256 0 1 2 1
A 770 7 0 14 1 10 1
A 777 7 258 0 1 2 1
A 776 7 0 14 1 10 1
A 783 7 260 0 1 2 1
A 782 7 0 14 1 10 1
A 789 7 262 0 1 2 1
A 788 7 0 14 1 10 1
A 795 7 264 0 1 2 1
A 794 7 0 14 1 10 1
A 801 7 266 0 1 2 1
A 800 7 0 14 1 10 1
A 809 7 268 0 1 2 1
A 808 7 0 164 1 10 1
A 817 7 270 0 1 2 1
A 816 7 0 164 1 10 1
A 825 7 272 0 1 2 1
A 824 7 0 164 1 10 0
T 1069 322 0 0 0 0
A 1077 7 334 0 1 2 1
A 1076 7 0 367 1 10 0
Z
