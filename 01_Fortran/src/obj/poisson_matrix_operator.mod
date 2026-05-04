V34 :0x24 poisson_matrix_operator
27 poisson_matrix_operator.f90 S624 0
05/04/2026  18:30:15
use geometry private
use matrix private
enduse
D 58 23 10 3 12 24 1 1 0 0 1
 10 13 11 14 13 15
 10 16 17 18 16 19
 10 20 21 22 20 23
D 61 23 10 3 25 37 1 1 0 0 1
 10 26 11 27 26 28
 10 29 30 31 29 32
 10 33 34 35 33 36
D 64 23 18 1 10 39 0 0 0 0 0
 10 38 11 10 38 39
D 91 26 708 1448 707 7
D 154 22 7
D 156 22 7
D 158 22 7
D 160 22 7
D 162 22 7
D 164 22 7
D 166 22 7
D 168 22 7
D 170 22 7
D 175 26 775 2288 774 7
D 259 22 7
D 261 22 7
D 263 22 7
D 265 22 7
D 267 22 7
D 269 22 7
D 271 22 7
D 273 22 7
D 275 22 7
D 277 22 7
D 279 22 7
D 281 22 7
D 331 26 1104 304 1103 7
D 343 22 7
D 408 26 775 2288 774 7
D 414 26 1104 304 1103 7
D 420 23 10 3 419 431 1 1 0 0 1
 10 420 11 421 420 422
 10 423 424 425 423 426
 10 427 428 429 427 430
D 423 23 10 3 432 444 1 1 0 0 1
 10 433 11 434 433 435
 10 436 437 438 436 439
 10 440 441 442 440 443
D 426 23 18 1 10 39 0 0 0 0 0
 10 38 11 10 38 39
S 624 24 0 0 0 9 1 0 5013 10005 0 A 0 0 0 0 B 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 poisson_matrix_operator
S 625 23 5 0 0 0 633 624 5037 0 0 A 0 0 0 0 B 0 7 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 vv_dot_3d_matrix
S 626 1 3 2 0 10 1 625 5054 4 3000 A 0 0 0 0 B 0 7 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 result
S 627 7 3 1 0 58 1 625 5061 20000004 10003000 A 0 0 0 0 B 0 7 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 x
S 628 7 3 1 0 61 1 625 5063 20000004 10003000 A 0 0 0 0 B 0 7 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 y
S 629 1 3 1 0 6 1 625 5065 4 3000 A 0 0 0 0 B 0 7 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 nx
S 630 1 3 1 0 6 1 625 5068 4 3000 A 0 0 0 0 B 0 7 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ny
S 631 1 3 1 0 6 1 625 5071 4 3000 A 0 0 0 0 B 0 7 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 nz
S 632 7 3 1 0 64 1 625 5074 800004 3000 A 0 0 0 0 B 0 7 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 is_serial
S 633 14 5 0 0 0 1 625 5037 20000000 400000 A 0 0 0 0 B 0 7 0 0 0 0 0 2 7 0 0 0 0 0 0 0 0 0 0 0 0 7 0 624 0 0 0 0 vv_dot_3d_matrix vv_dot_3d_matrix 
F 633 7 626 627 628 629 630 631 632
S 634 3 0 0 0 7 1 1 0 0 0 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 2 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 7
S 635 3 0 0 0 7 1 1 0 0 0 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 3 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 7
S 636 6 1 0 0 7 1 625 5084 40800006 3000 A 0 0 0 0 B 0 15 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 z_b_0
S 637 6 1 0 0 7 1 625 5090 40800006 3000 A 0 0 0 0 B 0 15 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 z_b_1
S 638 6 1 0 0 7 1 625 5096 40800006 3000 A 0 0 0 0 B 0 15 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 z_b_2
S 639 6 1 0 0 7 1 625 5102 40800006 3000 A 0 0 0 0 B 0 15 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 z_b_3
S 640 6 1 0 0 7 1 625 5108 40800006 3000 A 0 0 0 0 B 0 15 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 z_b_4
S 641 6 1 0 0 7 1 625 5114 40800006 3000 A 0 0 0 0 B 0 15 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 z_b_5
S 642 6 1 0 0 7 1 625 5120 40800006 3000 A 0 0 0 0 B 0 15 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 z_b_6
S 643 6 1 0 0 7 1 625 5126 40800006 3000 A 0 0 0 0 B 0 15 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 z_b_7
S 644 6 1 0 0 7 1 625 5132 40800006 3000 A 0 0 0 0 B 0 15 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 z_b_8
S 645 6 1 0 0 7 1 625 5138 40800006 3000 A 0 0 0 0 B 0 15 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 z_b_9
S 646 6 1 0 0 7 1 625 5144 40800006 3000 A 0 0 0 0 B 0 15 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 z_e_1320
S 647 6 1 0 0 7 1 625 5153 40800006 3000 A 0 0 0 0 B 0 15 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 z_e_1323
S 648 6 1 0 0 7 1 625 5162 40800006 3000 A 0 0 0 0 B 0 15 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 z_e_1326
S 649 6 1 0 0 7 1 625 5171 40800006 3000 A 0 0 0 0 B 0 15 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 z_b_10
S 650 6 1 0 0 7 1 625 5178 40800006 3000 A 0 0 0 0 B 0 15 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 z_b_11
S 651 6 1 0 0 7 1 625 5185 40800006 3000 A 0 0 0 0 B 0 15 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 z_b_12
S 652 6 1 0 0 7 1 625 5192 40800006 3000 A 0 0 0 0 B 0 15 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 z_b_13
S 653 6 1 0 0 7 1 625 5199 40800006 3000 A 0 0 0 0 B 0 15 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 z_b_14
S 654 6 1 0 0 7 1 625 5206 40800006 3000 A 0 0 0 0 B 0 15 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 z_b_15
S 655 6 1 0 0 7 1 625 5213 40800006 3000 A 0 0 0 0 B 0 15 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 z_b_16
S 656 6 1 0 0 7 1 625 5220 40800006 3000 A 0 0 0 0 B 0 15 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 z_b_17
S 657 6 1 0 0 7 1 625 5227 40800006 3000 A 0 0 0 0 B 0 15 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 z_b_18
S 658 6 1 0 0 7 1 625 5234 40800006 3000 A 0 0 0 0 B 0 15 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 z_b_19
S 659 6 1 0 0 7 1 625 5241 40800006 3000 A 0 0 0 0 B 0 15 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 z_e_1339
S 660 6 1 0 0 7 1 625 5250 40800006 3000 A 0 0 0 0 B 0 15 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 z_e_1342
S 661 6 1 0 0 7 1 625 5259 40800006 3000 A 0 0 0 0 B 0 15 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 z_e_1345
S 695 3 0 0 0 7 1 1 0 0 0 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 17 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 7
S 701 3 0 0 0 7 1 1 0 0 0 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 29 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 7
R 707 25 1 geometry domain
R 708 5 2 geometry nx domain
R 709 5 3 geometry ny domain
R 710 5 4 geometry nz domain
R 711 5 5 geometry lx domain
R 712 5 6 geometry ly domain
R 713 5 7 geometry lz domain
R 714 5 8 geometry ox domain
R 715 5 9 geometry oy domain
R 716 5 10 geometry oz domain
R 718 5 12 geometry dxm domain
R 719 5 13 geometry dxm$sd domain
R 720 5 14 geometry dxm$p domain
R 721 5 15 geometry dxm$o domain
R 723 5 17 geometry dym domain
R 725 5 19 geometry dym$sd domain
R 726 5 20 geometry dym$p domain
R 727 5 21 geometry dym$o domain
R 729 5 23 geometry dzm domain
R 731 5 25 geometry dzm$sd domain
R 732 5 26 geometry dzm$p domain
R 733 5 27 geometry dzm$o domain
R 736 5 30 geometry dxg domain
R 737 5 31 geometry dxg$sd domain
R 738 5 32 geometry dxg$p domain
R 739 5 33 geometry dxg$o domain
R 741 5 35 geometry dyg domain
R 743 5 37 geometry dyg$sd domain
R 744 5 38 geometry dyg$p domain
R 745 5 39 geometry dyg$o domain
R 747 5 41 geometry dzg domain
R 749 5 43 geometry dzg$sd domain
R 750 5 44 geometry dzg$p domain
R 751 5 45 geometry dzg$o domain
R 754 5 48 geometry xg domain
R 755 5 49 geometry xg$sd domain
R 756 5 50 geometry xg$p domain
R 757 5 51 geometry xg$o domain
R 759 5 53 geometry yg domain
R 761 5 55 geometry yg$sd domain
R 762 5 56 geometry yg$p domain
R 763 5 57 geometry yg$o domain
R 765 5 59 geometry zg domain
R 767 5 61 geometry zg$sd domain
R 768 5 62 geometry zg$p domain
R 769 5 63 geometry zg$o domain
R 771 5 65 geometry is_periodic domain
R 774 25 68 geometry subdomain
R 775 5 69 geometry nx subdomain
R 776 5 70 geometry ny subdomain
R 777 5 71 geometry nz subdomain
R 778 5 72 geometry lx subdomain
R 779 5 73 geometry ly subdomain
R 780 5 74 geometry lz subdomain
R 781 5 75 geometry ox subdomain
R 782 5 76 geometry oy subdomain
R 783 5 77 geometry oz subdomain
R 785 5 79 geometry dxm subdomain
R 786 5 80 geometry dxm$sd subdomain
R 787 5 81 geometry dxm$p subdomain
R 788 5 82 geometry dxm$o subdomain
R 791 5 85 geometry dym subdomain
R 792 5 86 geometry dym$sd subdomain
R 793 5 87 geometry dym$p subdomain
R 794 5 88 geometry dym$o subdomain
R 797 5 91 geometry dzm subdomain
R 798 5 92 geometry dzm$sd subdomain
R 799 5 93 geometry dzm$p subdomain
R 800 5 94 geometry dzm$o subdomain
R 803 5 97 geometry dxg subdomain
R 804 5 98 geometry dxg$sd subdomain
R 805 5 99 geometry dxg$p subdomain
R 806 5 100 geometry dxg$o subdomain
R 809 5 103 geometry dyg subdomain
R 810 5 104 geometry dyg$sd subdomain
R 811 5 105 geometry dyg$p subdomain
R 812 5 106 geometry dyg$o subdomain
R 815 5 109 geometry dzg subdomain
R 816 5 110 geometry dzg$sd subdomain
R 817 5 111 geometry dzg$p subdomain
R 818 5 112 geometry dzg$o subdomain
R 821 5 115 geometry xg subdomain
R 822 5 116 geometry xg$sd subdomain
R 823 5 117 geometry xg$p subdomain
R 824 5 118 geometry xg$o subdomain
R 827 5 121 geometry yg subdomain
R 828 5 122 geometry yg$sd subdomain
R 829 5 123 geometry yg$p subdomain
R 830 5 124 geometry yg$o subdomain
R 833 5 127 geometry zg subdomain
R 834 5 128 geometry zg$sd subdomain
R 835 5 129 geometry zg$p subdomain
R 836 5 130 geometry zg$o subdomain
R 841 5 135 geometry x subdomain
R 842 5 136 geometry x$sd subdomain
R 843 5 137 geometry x$p subdomain
R 844 5 138 geometry x$o subdomain
R 846 5 140 geometry b subdomain
R 850 5 144 geometry b$sd subdomain
R 851 5 145 geometry b$p subdomain
R 852 5 146 geometry b$o subdomain
R 854 5 148 geometry r subdomain
R 858 5 152 geometry r$sd subdomain
R 859 5 153 geometry r$p subdomain
R 860 5 154 geometry r$o subdomain
R 862 5 156 geometry is_periodic subdomain
R 863 5 157 geometry is_aggregated subdomain
R 864 5 158 geometry ista subdomain
R 865 5 159 geometry iend subdomain
R 866 5 160 geometry jsta subdomain
R 867 5 161 geometry jend subdomain
R 868 5 162 geometry ksta subdomain
R 869 5 163 geometry kend subdomain
R 870 5 164 geometry ddt_yz_plane_x0 subdomain
R 871 5 165 geometry ddt_yz_plane_x1 subdomain
R 872 5 166 geometry ddt_yz_plane_xn subdomain
R 873 5 167 geometry ddt_yz_plane_xn1 subdomain
R 874 5 168 geometry ddt_xz_plane_y0 subdomain
R 875 5 169 geometry ddt_xz_plane_y1 subdomain
R 876 5 170 geometry ddt_xz_plane_yn subdomain
R 877 5 171 geometry ddt_xz_plane_yn1 subdomain
R 878 5 172 geometry ddt_xy_plane_z0 subdomain
R 879 5 173 geometry ddt_xy_plane_z1 subdomain
R 880 5 174 geometry ddt_xy_plane_zn subdomain
R 881 5 175 geometry ddt_xy_plane_zn1 subdomain
R 882 5 176 geometry ddt_inner_domain subdomain
R 883 5 177 geometry is_x0_boundary subdomain
R 884 5 178 geometry is_x1_boundary subdomain
R 885 5 179 geometry is_y0_boundary subdomain
R 886 5 180 geometry is_y1_boundary subdomain
R 887 5 181 geometry is_z0_boundary subdomain
R 888 5 182 geometry is_z1_boundary subdomain
S 1100 3 0 0 0 7 1 1 0 0 0 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 35 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 7
R 1103 25 1 matrix matrix_poisson
R 1104 5 2 matrix dof matrix_poisson
R 1105 5 3 matrix coeff matrix_poisson
R 1110 5 8 matrix coeff$sd matrix_poisson
R 1111 5 9 matrix coeff$p matrix_poisson
R 1112 5 10 matrix coeff$o matrix_poisson
S 1125 23 5 0 0 0 1131 624 7715 0 0 A 0 0 0 0 B 0 57 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 mv_mul_poisson_matrix
S 1126 7 3 2 0 420 1 1125 5063 20000004 10003000 A 0 0 0 0 B 0 57 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 y
S 1127 1 3 1 0 414 1 1125 7682 4 3000 A 0 0 0 0 B 0 57 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 a_poisson
S 1128 7 3 3 0 423 1 1125 5061 20000004 10003000 A 0 0 0 0 B 0 57 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 x
S 1129 1 3 1 0 408 1 1125 7737 4 3000 A 0 0 0 0 B 0 57 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 dm
S 1130 7 3 1 0 426 1 1125 5074 800004 3000 A 0 0 0 0 B 0 57 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 is_serial
S 1131 14 5 0 0 0 1 1125 7715 20000000 400000 A 0 0 0 0 B 0 57 0 0 0 0 0 69 5 0 0 0 0 0 0 0 0 0 0 0 0 57 0 624 0 0 0 0 mv_mul_poisson_matrix mv_mul_poisson_matrix 
F 1131 5 1126 1127 1128 1129 1130
S 1132 6 1 0 0 7 1 1125 6778 40800006 3000 A 0 0 0 0 B 0 65 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 z_b_0_1
S 1133 6 1 0 0 7 1 1125 6786 40800006 3000 A 0 0 0 0 B 0 65 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 z_b_1_1
S 1134 6 1 0 0 7 1 1125 6794 40800006 3000 A 0 0 0 0 B 0 65 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 z_b_2_1
S 1135 6 1 0 0 7 1 1125 6802 40800006 3000 A 0 0 0 0 B 0 65 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 z_b_3_1
S 1136 6 1 0 0 7 1 1125 6810 40800006 3000 A 0 0 0 0 B 0 65 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 z_b_4_1
S 1137 6 1 0 0 7 1 1125 6818 40800006 3000 A 0 0 0 0 B 0 65 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 z_b_5_1
S 1138 6 1 0 0 7 1 1125 6826 40800006 3000 A 0 0 0 0 B 0 65 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 z_b_6_1
S 1139 6 1 0 0 7 1 1125 6834 40800006 3000 A 0 0 0 0 B 0 65 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 z_b_7_1
S 1140 6 1 0 0 7 1 1125 6842 40800006 3000 A 0 0 0 0 B 0 65 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 z_b_8_1
S 1141 6 1 0 0 7 1 1125 6850 40800006 3000 A 0 0 0 0 B 0 65 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 z_b_9_1
S 1142 6 1 0 0 7 1 1125 7740 40800006 3000 A 0 0 0 0 B 0 65 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 z_e_1724
S 1143 6 1 0 0 7 1 1125 7749 40800006 3000 A 0 0 0 0 B 0 65 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 z_e_1727
S 1144 6 1 0 0 7 1 1125 7758 40800006 3000 A 0 0 0 0 B 0 65 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 z_e_1730
S 1145 6 1 0 0 7 1 1125 7085 40800006 3000 A 0 0 0 0 B 0 67 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 z_b_10_1
S 1146 6 1 0 0 7 1 1125 7094 40800006 3000 A 0 0 0 0 B 0 67 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 z_b_11_1
S 1147 6 1 0 0 7 1 1125 7103 40800006 3000 A 0 0 0 0 B 0 67 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 z_b_12_1
S 1148 6 1 0 0 7 1 1125 7112 40800006 3000 A 0 0 0 0 B 0 67 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 z_b_13_1
S 1149 6 1 0 0 7 1 1125 7121 40800006 3000 A 0 0 0 0 B 0 67 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 z_b_14_1
S 1150 6 1 0 0 7 1 1125 7130 40800006 3000 A 0 0 0 0 B 0 67 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 z_b_15_1
S 1151 6 1 0 0 7 1 1125 7139 40800006 3000 A 0 0 0 0 B 0 67 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 z_b_16_1
S 1152 6 1 0 0 7 1 1125 7148 40800006 3000 A 0 0 0 0 B 0 67 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 z_b_17_1
S 1153 6 1 0 0 7 1 1125 7157 40800006 3000 A 0 0 0 0 B 0 67 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 z_b_18_1
S 1154 6 1 0 0 7 1 1125 7166 40800006 3000 A 0 0 0 0 B 0 67 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 z_b_19_1
S 1155 6 1 0 0 7 1 1125 7767 40800006 3000 A 0 0 0 0 B 0 67 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 z_e_1743
S 1156 6 1 0 0 7 1 1125 7776 40800006 3000 A 0 0 0 0 B 0 67 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 z_e_1746
S 1157 6 1 0 0 7 1 1125 7785 40800006 3000 A 0 0 0 0 B 0 67 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 z_e_1749
A 12 1 0 0 0 7 645 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
A 13 1 0 0 0 7 637 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
A 14 1 0 0 0 7 636 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
A 15 1 0 0 0 7 646 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
A 16 1 0 0 0 7 640 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
A 17 1 0 0 0 7 638 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
A 18 1 0 0 0 7 639 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
A 19 1 0 0 0 7 647 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
A 20 1 0 0 0 7 643 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
A 21 1 0 0 0 7 641 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
A 22 1 0 0 0 7 642 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
A 23 1 0 0 0 7 648 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
A 24 1 0 0 0 7 644 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
A 25 1 0 0 0 7 658 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
A 26 1 0 0 0 7 650 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
A 27 1 0 0 0 7 649 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
A 28 1 0 0 0 7 659 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
A 29 1 0 0 0 7 653 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
A 30 1 0 0 0 7 651 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
A 31 1 0 0 0 7 652 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
A 32 1 0 0 0 7 660 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
A 33 1 0 0 0 7 656 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
A 34 1 0 0 0 7 654 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
A 35 1 0 0 0 7 655 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
A 36 1 0 0 0 7 661 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
A 37 1 0 0 0 7 657 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
A 38 2 0 0 0 7 634 0 0 0 38 0 0 0 0 0 0 0 0 0 0 0
A 39 2 0 0 0 7 635 0 0 0 39 0 0 0 0 0 0 0 0 0 0 0
A 40 2 0 0 0 7 695 0 0 0 40 0 0 0 0 0 0 0 0 0 0 0
A 190 2 0 0 0 7 701 0 0 0 190 0 0 0 0 0 0 0 0 0 0 0
A 393 2 0 0 0 7 1100 0 0 0 393 0 0 0 0 0 0 0 0 0 0 0
A 419 1 0 0 0 7 1141 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
A 420 1 0 0 0 7 1133 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
A 421 1 0 0 0 7 1132 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
A 422 1 0 0 0 7 1142 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
A 423 1 0 0 0 7 1136 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
A 424 1 0 0 0 7 1134 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
A 425 1 0 0 0 7 1135 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
A 426 1 0 0 0 7 1143 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
A 427 1 0 0 0 7 1139 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
A 428 1 0 0 0 7 1137 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
A 429 1 0 0 0 7 1138 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
A 430 1 0 0 0 7 1144 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
A 431 1 0 0 0 7 1140 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
A 432 1 0 0 0 7 1154 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
A 433 1 0 0 0 7 1146 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
A 434 1 0 0 0 7 1145 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
A 435 1 0 0 0 7 1155 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
A 436 1 0 0 196 7 1149 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
A 437 1 0 0 0 7 1147 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
A 438 1 0 0 0 7 1148 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
A 439 1 0 0 0 7 1156 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
A 440 1 0 0 0 7 1152 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
A 441 1 0 0 0 7 1150 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
A 442 1 0 0 0 7 1151 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
A 443 1 0 0 0 7 1157 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
A 444 1 0 0 0 7 1153 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
Z
T 707 91 0 0 0 0
A 720 7 154 0 1 2 1
A 719 7 0 40 1 10 1
A 726 7 156 0 1 2 1
A 725 7 0 40 1 10 1
A 732 7 158 0 1 2 1
A 731 7 0 40 1 10 1
A 738 7 160 0 1 2 1
A 737 7 0 40 1 10 1
A 744 7 162 0 1 2 1
A 743 7 0 40 1 10 1
A 750 7 164 0 1 2 1
A 749 7 0 40 1 10 1
A 756 7 166 0 1 2 1
A 755 7 0 40 1 10 1
A 762 7 168 0 1 2 1
A 761 7 0 40 1 10 1
A 768 7 170 0 1 2 1
A 767 7 0 40 1 10 0
T 774 175 0 0 0 0
A 787 7 259 0 1 2 1
A 786 7 0 40 1 10 1
A 793 7 261 0 1 2 1
A 792 7 0 40 1 10 1
A 799 7 263 0 1 2 1
A 798 7 0 40 1 10 1
A 805 7 265 0 1 2 1
A 804 7 0 40 1 10 1
A 811 7 267 0 1 2 1
A 810 7 0 40 1 10 1
A 817 7 269 0 1 2 1
A 816 7 0 40 1 10 1
A 823 7 271 0 1 2 1
A 822 7 0 40 1 10 1
A 829 7 273 0 1 2 1
A 828 7 0 40 1 10 1
A 835 7 275 0 1 2 1
A 834 7 0 40 1 10 1
A 843 7 277 0 1 2 1
A 842 7 0 190 1 10 1
A 851 7 279 0 1 2 1
A 850 7 0 190 1 10 1
A 859 7 281 0 1 2 1
A 858 7 0 190 1 10 0
T 1103 331 0 0 0 0
A 1111 7 343 0 1 2 1
A 1110 7 0 393 1 10 0
Z
