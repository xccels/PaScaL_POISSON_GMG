V34 :0x24 matrix
10 matrix.f90 S624 0
05/04/2026  18:30:15
use geometry private
enduse
D 58 26 627 304 625 7
D 64 23 10 4 55 53 0 1 0 0 1
 21 25 45 21 25 23
 27 31 47 27 31 29
 33 37 49 33 37 35
 39 43 51 39 43 41
D 67 23 7 1 0 18 0 0 0 0 0
 0 18 0 11 18 0
D 70 22 7
D 72 23 7 1 0 11 0 0 0 0 0
 0 11 0 11 11 0
D 99 26 694 1448 693 7
D 162 22 7
D 164 22 7
D 166 22 7
D 168 22 7
D 170 22 7
D 172 22 7
D 174 22 7
D 176 22 7
D 178 22 7
D 183 26 761 2288 760 7
D 267 22 7
D 269 22 7
D 271 22 7
D 273 22 7
D 275 22 7
D 277 22 7
D 279 22 7
D 281 22 7
D 283 22 7
D 285 22 7
D 287 22 7
D 289 22 7
D 339 26 761 2288 760 7
S 624 24 0 0 0 6 1 0 5013 10005 0 A 0 0 0 0 B 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 matrix
S 625 25 0 0 0 58 1 624 5020 1000000c 800050 A 0 0 0 0 B 0 5 0 0 0 0 0 0 0 0 0 655 0 0 0 0 0 0 0 654 0 0 0 624 0 0 0 0 matrix_poisson
S 627 5 0 0 0 6 629 624 5035 800004 0 A 0 0 0 0 B 0 6 0 0 0 0 0 0 58 0 0 0 0 0 0 0 0 0 0 0 1 627 0 624 0 0 0 0 dof
S 629 5 6 0 0 64 636 624 5039 10a00004 51 A 0 0 0 0 B 0 7 0 0 0 8 636 0 58 0 638 0 0 0 0 0 0 0 0 635 627 629 637 624 0 0 0 0 coeff
S 630 6 4 0 0 7 631 624 5045 40800006 0 A 0 0 0 0 B 0 7 0 0 0 0 0 0 0 0 0 0 656 0 0 0 0 0 0 0 0 0 0 624 0 0 0 0 z_b_0
S 631 6 4 0 0 7 632 624 5051 40800006 0 A 0 0 0 0 B 0 7 0 0 0 8 0 0 0 0 0 0 656 0 0 0 0 0 0 0 0 0 0 624 0 0 0 0 z_b_1
S 632 6 4 0 0 7 633 624 5057 40800006 0 A 0 0 0 0 B 0 7 0 0 0 16 0 0 0 0 0 0 656 0 0 0 0 0 0 0 0 0 0 624 0 0 0 0 z_b_2
S 633 6 4 0 0 7 1 624 5063 40800006 0 A 0 0 0 0 B 0 7 0 0 0 24 0 0 0 0 0 0 656 0 0 0 0 0 0 0 0 0 0 624 0 0 0 0 z_b_3
S 634 3 0 0 0 7 1 1 0 0 0 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 35 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 7
S 635 5 1 0 0 67 1 624 5069 40822004 1020 A 0 0 0 0 B 0 7 0 0 0 24 0 0 58 0 0 0 0 0 0 0 0 0 0 0 637 635 0 624 0 0 0 0 coeff$sd
S 636 5 0 0 0 7 637 624 5078 40802001 1020 A 0 0 0 0 B 0 7 0 0 0 8 0 0 58 0 0 0 0 0 0 0 0 0 0 0 629 636 0 624 0 0 0 0 coeff$p
S 637 5 0 0 0 7 635 624 5086 40802000 1020 A 0 0 0 0 B 0 7 0 0 0 16 0 0 58 0 0 0 0 0 0 0 0 0 0 0 636 637 0 624 0 0 0 0 coeff$o
S 638 22 1 0 0 9 1 624 5094 40000000 1000 A 0 0 0 0 B 0 7 0 0 0 0 0 629 0 0 0 0 635 0 0 0 0 0 0 0 0 0 0 624 0 0 0 0 coeff$arrdsc
S 639 3 0 0 0 7 1 1 0 0 0 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 11 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 7
S 640 3 0 0 0 7 1 1 0 0 0 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 12 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 7
S 642 3 0 0 0 7 1 1 0 0 0 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 17 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 7
S 643 3 0 0 0 7 1 1 0 0 0 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 18 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 7
S 644 3 0 0 0 7 1 1 0 0 0 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 23 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 7
S 645 3 0 0 0 7 1 1 0 0 0 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 24 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 7
S 646 3 0 0 0 7 1 1 0 0 0 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 29 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 7
S 647 3 0 0 0 7 1 1 0 0 0 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 30 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 7
S 648 3 0 0 0 7 1 1 0 0 0 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 15 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 7
S 649 3 0 0 0 7 1 1 0 0 0 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 21 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 7
S 650 3 0 0 0 7 1 1 0 0 0 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 27 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 7
S 651 3 0 0 0 7 1 1 0 0 0 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 33 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 7
S 652 3 0 0 0 7 1 1 0 0 0 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 7 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 7
S 653 3 0 0 0 7 1 1 0 0 0 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 8 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 7
S 654 8 5 0 0 72 1 624 5107 40822004 1220 A 0 0 0 0 B 0 8 0 0 0 0 0 58 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 624 0 0 0 0 matrix$matrix_poisson$$td
S 655 6 4 0 0 58 1 624 5133 80004e 0 A 0 0 0 0 B 800 8 0 0 0 0 0 0 0 0 0 0 657 0 0 0 0 0 0 0 0 0 0 624 0 0 0 0 ._dtInit0058
S 656 11 0 0 0 9 1 624 5146 40800000 805000 A 0 0 0 0 B 0 10 0 0 0 32 0 0 630 633 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 _matrix$2
S 657 11 0 0 0 9 656 624 5156 40800000 805000 A 0 0 0 0 B 0 10 0 0 0 304 0 0 655 655 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 _matrix$8
R 693 25 1 geometry domain
R 694 5 2 geometry nx domain
R 695 5 3 geometry ny domain
R 696 5 4 geometry nz domain
R 697 5 5 geometry lx domain
R 698 5 6 geometry ly domain
R 699 5 7 geometry lz domain
R 700 5 8 geometry ox domain
R 701 5 9 geometry oy domain
R 702 5 10 geometry oz domain
R 704 5 12 geometry dxm domain
R 705 5 13 geometry dxm$sd domain
R 706 5 14 geometry dxm$p domain
R 707 5 15 geometry dxm$o domain
R 709 5 17 geometry dym domain
R 711 5 19 geometry dym$sd domain
R 712 5 20 geometry dym$p domain
R 713 5 21 geometry dym$o domain
R 715 5 23 geometry dzm domain
R 717 5 25 geometry dzm$sd domain
R 718 5 26 geometry dzm$p domain
R 719 5 27 geometry dzm$o domain
R 722 5 30 geometry dxg domain
R 723 5 31 geometry dxg$sd domain
R 724 5 32 geometry dxg$p domain
R 725 5 33 geometry dxg$o domain
R 727 5 35 geometry dyg domain
R 729 5 37 geometry dyg$sd domain
R 730 5 38 geometry dyg$p domain
R 731 5 39 geometry dyg$o domain
R 733 5 41 geometry dzg domain
R 735 5 43 geometry dzg$sd domain
R 736 5 44 geometry dzg$p domain
R 737 5 45 geometry dzg$o domain
R 740 5 48 geometry xg domain
R 741 5 49 geometry xg$sd domain
R 742 5 50 geometry xg$p domain
R 743 5 51 geometry xg$o domain
R 745 5 53 geometry yg domain
R 747 5 55 geometry yg$sd domain
R 748 5 56 geometry yg$p domain
R 749 5 57 geometry yg$o domain
R 751 5 59 geometry zg domain
R 753 5 61 geometry zg$sd domain
R 754 5 62 geometry zg$p domain
R 755 5 63 geometry zg$o domain
R 757 5 65 geometry is_periodic domain
R 760 25 68 geometry subdomain
R 761 5 69 geometry nx subdomain
R 762 5 70 geometry ny subdomain
R 763 5 71 geometry nz subdomain
R 764 5 72 geometry lx subdomain
R 765 5 73 geometry ly subdomain
R 766 5 74 geometry lz subdomain
R 767 5 75 geometry ox subdomain
R 768 5 76 geometry oy subdomain
R 769 5 77 geometry oz subdomain
R 771 5 79 geometry dxm subdomain
R 772 5 80 geometry dxm$sd subdomain
R 773 5 81 geometry dxm$p subdomain
R 774 5 82 geometry dxm$o subdomain
R 777 5 85 geometry dym subdomain
R 778 5 86 geometry dym$sd subdomain
R 779 5 87 geometry dym$p subdomain
R 780 5 88 geometry dym$o subdomain
R 783 5 91 geometry dzm subdomain
R 784 5 92 geometry dzm$sd subdomain
R 785 5 93 geometry dzm$p subdomain
R 786 5 94 geometry dzm$o subdomain
R 789 5 97 geometry dxg subdomain
R 790 5 98 geometry dxg$sd subdomain
R 791 5 99 geometry dxg$p subdomain
R 792 5 100 geometry dxg$o subdomain
R 795 5 103 geometry dyg subdomain
R 796 5 104 geometry dyg$sd subdomain
R 797 5 105 geometry dyg$p subdomain
R 798 5 106 geometry dyg$o subdomain
R 801 5 109 geometry dzg subdomain
R 802 5 110 geometry dzg$sd subdomain
R 803 5 111 geometry dzg$p subdomain
R 804 5 112 geometry dzg$o subdomain
R 807 5 115 geometry xg subdomain
R 808 5 116 geometry xg$sd subdomain
R 809 5 117 geometry xg$p subdomain
R 810 5 118 geometry xg$o subdomain
R 813 5 121 geometry yg subdomain
R 814 5 122 geometry yg$sd subdomain
R 815 5 123 geometry yg$p subdomain
R 816 5 124 geometry yg$o subdomain
R 819 5 127 geometry zg subdomain
R 820 5 128 geometry zg$sd subdomain
R 821 5 129 geometry zg$p subdomain
R 822 5 130 geometry zg$o subdomain
R 827 5 135 geometry x subdomain
R 828 5 136 geometry x$sd subdomain
R 829 5 137 geometry x$p subdomain
R 830 5 138 geometry x$o subdomain
R 832 5 140 geometry b subdomain
R 836 5 144 geometry b$sd subdomain
R 837 5 145 geometry b$p subdomain
R 838 5 146 geometry b$o subdomain
R 840 5 148 geometry r subdomain
R 844 5 152 geometry r$sd subdomain
R 845 5 153 geometry r$p subdomain
R 846 5 154 geometry r$o subdomain
R 848 5 156 geometry is_periodic subdomain
R 849 5 157 geometry is_aggregated subdomain
R 850 5 158 geometry ista subdomain
R 851 5 159 geometry iend subdomain
R 852 5 160 geometry jsta subdomain
R 853 5 161 geometry jend subdomain
R 854 5 162 geometry ksta subdomain
R 855 5 163 geometry kend subdomain
R 856 5 164 geometry ddt_yz_plane_x0 subdomain
R 857 5 165 geometry ddt_yz_plane_x1 subdomain
R 858 5 166 geometry ddt_yz_plane_xn subdomain
R 859 5 167 geometry ddt_yz_plane_xn1 subdomain
R 860 5 168 geometry ddt_xz_plane_y0 subdomain
R 861 5 169 geometry ddt_xz_plane_y1 subdomain
R 862 5 170 geometry ddt_xz_plane_yn subdomain
R 863 5 171 geometry ddt_xz_plane_yn1 subdomain
R 864 5 172 geometry ddt_xy_plane_z0 subdomain
R 865 5 173 geometry ddt_xy_plane_z1 subdomain
R 866 5 174 geometry ddt_xy_plane_zn subdomain
R 867 5 175 geometry ddt_xy_plane_zn1 subdomain
R 868 5 176 geometry ddt_inner_domain subdomain
R 869 5 177 geometry is_x0_boundary subdomain
R 870 5 178 geometry is_x1_boundary subdomain
R 871 5 179 geometry is_y0_boundary subdomain
R 872 5 180 geometry is_y1_boundary subdomain
R 873 5 181 geometry is_z0_boundary subdomain
R 874 5 182 geometry is_z1_boundary subdomain
S 1084 23 5 0 0 0 1087 624 7556 0 0 A 0 0 0 0 B 0 12 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 matrix_poisson_create
S 1085 1 3 3 0 58 1 1084 7578 4 3000 A 0 0 0 0 B 0 12 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 a_poisson
S 1086 1 3 1 0 339 1 1084 6661 4 3000 A 0 0 0 0 B 0 12 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 sdm
S 1087 14 5 0 0 0 1 1084 7556 0 400000 A 0 0 0 0 B 0 12 0 0 0 0 0 56 2 0 0 0 0 0 0 0 0 0 0 0 0 12 0 624 0 0 0 0 matrix_poisson_create matrix_poisson_create 
F 1087 2 1085 1086
S 1088 23 5 0 0 0 1090 624 7588 0 0 A 0 0 0 0 B 0 58 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 matrix_poisson_destroy
S 1089 1 3 3 0 58 1 1088 7578 4 3000 A 0 0 0 0 B 0 58 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 a_poisson
S 1090 14 5 0 0 0 1 1088 7588 0 400000 A 0 0 0 0 B 0 58 0 0 0 0 0 59 1 0 0 0 0 0 0 0 0 0 0 0 0 58 0 624 0 0 0 0 matrix_poisson_destroy matrix_poisson_destroy 
F 1090 1 1089
A 18 2 0 0 0 7 634 0 0 0 18 0 0 0 0 0 0 0 0 0 0 0
A 19 2 0 0 0 7 639 0 0 0 19 0 0 0 0 0 0 0 0 0 0 0
A 20 1 0 1 0 67 635 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
A 21 10 0 0 0 7 20 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0
X 1 19
A 22 2 0 0 0 7 640 0 0 0 22 0 0 0 0 0 0 0 0 0 0 0
A 23 10 0 0 21 7 20 4 0 0 0 0 0 0 0 0 0 0 0 0 0 0
X 1 22
A 24 4 0 0 0 7 23 0 11 0 0 0 0 2 0 0 0 0 0 0 0 0
A 25 4 0 0 0 7 21 0 24 0 0 0 0 1 0 0 0 0 0 0 0 0
A 26 2 0 0 0 7 642 0 0 0 26 0 0 0 0 0 0 0 0 0 0 0
A 27 10 0 0 23 7 20 7 0 0 0 0 0 0 0 0 0 0 0 0 0 0
X 1 26
A 28 2 0 0 0 7 643 0 0 0 28 0 0 0 0 0 0 0 0 0 0 0
A 29 10 0 0 27 7 20 10 0 0 0 0 0 0 0 0 0 0 0 0 0 0
X 1 28
A 30 4 0 0 0 7 29 0 11 0 0 0 0 2 0 0 0 0 0 0 0 0
A 31 4 0 0 0 7 27 0 30 0 0 0 0 1 0 0 0 0 0 0 0 0
A 32 2 0 0 0 7 644 0 0 0 32 0 0 0 0 0 0 0 0 0 0 0
A 33 10 0 0 29 7 20 13 0 0 0 0 0 0 0 0 0 0 0 0 0 0
X 1 32
A 34 2 0 0 0 7 645 0 0 0 34 0 0 0 0 0 0 0 0 0 0 0
A 35 10 0 0 33 7 20 16 0 0 0 0 0 0 0 0 0 0 0 0 0 0
X 1 34
A 36 4 0 0 0 7 35 0 11 0 0 0 0 2 0 0 0 0 0 0 0 0
A 37 4 0 0 0 7 33 0 36 0 0 0 0 1 0 0 0 0 0 0 0 0
A 38 2 0 0 0 7 646 0 0 0 38 0 0 0 0 0 0 0 0 0 0 0
A 39 10 0 0 35 7 20 19 0 0 0 0 0 0 0 0 0 0 0 0 0 0
X 1 38
A 40 2 0 0 0 7 647 0 0 0 40 0 0 0 0 0 0 0 0 0 0 0
A 41 10 0 0 39 7 20 22 0 0 0 0 0 0 0 0 0 0 0 0 0 0
X 1 40
A 42 4 0 0 0 7 41 0 11 0 0 0 0 2 0 0 0 0 0 0 0 0
A 43 4 0 0 0 7 39 0 42 0 0 0 0 1 0 0 0 0 0 0 0 0
A 44 2 0 0 0 7 648 0 0 0 44 0 0 0 0 0 0 0 0 0 0 0
A 45 10 0 0 41 7 20 25 0 0 0 0 0 0 0 0 0 0 0 0 0 0
X 1 44
A 46 2 0 0 0 7 649 0 0 0 46 0 0 0 0 0 0 0 0 0 0 0
A 47 10 0 0 45 7 20 28 0 0 0 0 0 0 0 0 0 0 0 0 0 0
X 1 46
A 48 2 0 0 0 7 650 0 0 0 48 0 0 0 0 0 0 0 0 0 0 0
A 49 10 0 0 47 7 20 31 0 0 0 0 0 0 0 0 0 0 0 0 0 0
X 1 48
A 50 2 0 0 0 7 651 0 0 0 50 0 0 0 0 0 0 0 0 0 0 0
A 51 10 0 0 49 7 20 34 0 0 0 0 0 0 0 0 0 0 0 0 0 0
X 1 50
A 52 2 0 0 0 7 652 0 0 0 52 0 0 0 0 0 0 0 0 0 0 0
A 53 10 0 0 51 7 20 37 0 0 0 0 0 0 0 0 0 0 0 0 0 0
X 1 52
A 54 2 0 0 0 7 653 0 0 0 54 0 0 0 0 0 0 0 0 0 0 0
A 55 10 0 0 53 7 20 40 0 0 0 0 0 0 0 0 0 0 0 0 0 0
X 1 54
Z
T 625 58 0 0 0 0
A 636 7 70 0 1 2 1
A 635 7 0 18 1 10 0
T 693 99 0 0 0 0
A 706 7 162 0 1 2 1
A 705 7 0 26 1 10 1
A 712 7 164 0 1 2 1
A 711 7 0 26 1 10 1
A 718 7 166 0 1 2 1
A 717 7 0 26 1 10 1
A 724 7 168 0 1 2 1
A 723 7 0 26 1 10 1
A 730 7 170 0 1 2 1
A 729 7 0 26 1 10 1
A 736 7 172 0 1 2 1
A 735 7 0 26 1 10 1
A 742 7 174 0 1 2 1
A 741 7 0 26 1 10 1
A 748 7 176 0 1 2 1
A 747 7 0 26 1 10 1
A 754 7 178 0 1 2 1
A 753 7 0 26 1 10 0
T 760 183 0 0 0 0
A 773 7 267 0 1 2 1
A 772 7 0 26 1 10 1
A 779 7 269 0 1 2 1
A 778 7 0 26 1 10 1
A 785 7 271 0 1 2 1
A 784 7 0 26 1 10 1
A 791 7 273 0 1 2 1
A 790 7 0 26 1 10 1
A 797 7 275 0 1 2 1
A 796 7 0 26 1 10 1
A 803 7 277 0 1 2 1
A 802 7 0 26 1 10 1
A 809 7 279 0 1 2 1
A 808 7 0 26 1 10 1
A 815 7 281 0 1 2 1
A 814 7 0 26 1 10 1
A 821 7 283 0 1 2 1
A 820 7 0 26 1 10 1
A 829 7 285 0 1 2 1
A 828 7 0 38 1 10 1
A 837 7 287 0 1 2 1
A 836 7 0 38 1 10 1
A 845 7 289 0 1 2 1
A 844 7 0 38 1 10 0
Z
