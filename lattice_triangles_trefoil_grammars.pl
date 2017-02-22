#!/usr/bin/perl -w
# HHC - 2016-06-16 - on men-pc1072 via men-pc4142
# HHC - 2016-06-17 - on men-pc2018 (ActivePerl v5.22.1)
# HHC - 2016-06-20 - on men-pc1072 via men-pc2018 (v5.10.1 from Debian)
# HHC - 2016-06-24 - on men-pc2018 Tk on Windows
# HHC - 2016-06-26 - at Orrington, Evanston near Chicago, IL
use strict;
use Tk;
use PDL;
use Set::Scalar;
use Memoize qw(memoize unmemoize);

my $pi = 4 * atan2(1,1);
my $epsilon1 = 1e-8;   # 0.01 um
my $epsilon2 = 1e-11;  # angle in radian within a 1km cube working space
my $epsilon3 = 1e-5,   # 0.01 mm
my $epsilon4 = 1e-16,  # pretty random for identical values

my (%coords, %arc, %weights, %curvature);
my %transformation;
my %shape;

### preliminaries and data


my ($p1_x, $p1_y,
           $p3_y,
           $c3_y,
    $c1_x, $c1_y,
    $centre_a2_x, $centre_a2_y) = &trefoil_control_polygon;

# $p1_x = x1 =  77.811_529_493_745_3
# $p1_y = y1 =  18.584_022_111_588_5
# $p3_y = y2 = 153.357_544_609_4
# $c3_y = c3 = 344.381_854_561_458
# $c1_x = x3 = 243.243_434_652_621
# $c1_y = y3 = -76.928_132_864_440_9
# a2_centre 55, 95.262_794_416_288_3  # not used, for reference only

$coords{p1} = transpose pdl [ $p1_x, $p1_y, 0];
$coords{c3} = transpose pdl [     0, $c3_y, 0];
$coords{p2} = transpose pdl [-$p1_x, $p1_y, 0];

#       p2 as well
$coords{c1} = transpose pdl [ $c1_x, $c1_y, 0];
$coords{p3} = transpose pdl [     0, $p3_y, 0];

#       p1 as well
$coords{c2} = transpose pdl [-$c1_x, $c1_y, 0];
#       p3 as well

$coords{m2} = transpose pdl [ 22.811_529_493_747_3, 76.678_772_304_708_3, 0];
$coords{m1} = transpose pdl [-22.811_529_493_747_3, 76.678_772_304_708_3, 0];
$coords{m3} = transpose pdl [  0                  , 37.168_044_223_172_5, 0];

&trefoil_rule_polygon( $p3_y - $coords{m3}->at(0,1) );

# U was called a3
# T was called a1
# V was called a2

$arc{U} = [qw/ p1 c3 p2 /];
$arc{T} = [qw/ p2 c1 p3 /];
$arc{V} = [qw/ p1 c2 p3 /];

while ( my $this = each %arc) {
    $weights{$this}   = [ &return_weights( @{$arc{$this}} ) ];
    $curvature{$this} = &return_curvature( @{$arc{$this}} );
    #print "arc $this : w1=$weights{$this}[1], 1/r=$curvature{$this}\n";
}

my @sets = (
    [ "U", "T", transpose pdl [ 20, 80, 0] ],
   #[ "T", "U", transpose pdl [ 20, 80, 0] ],
   #[ "U", "T", transpose pdl [-80, 20, 0] ],
   #[ "T", "U", transpose pdl [-80, 20, 0] ],

    [ "U", "V", transpose pdl [-20, 80, 0] ],
   #[ "V", "U", transpose pdl [-20, 80, 0] ],
   #[ "U", "V", transpose pdl [ 80, 20, 0] ],
   #[ "V", "U", transpose pdl [ 80, 20, 0] ],

    [ "T", "V", transpose pdl [  0, 40, 0] ],
   #[ "V", "T", transpose pdl [  0, 40, 0] ],
   #[ "T", "V", transpose pdl [  0,160, 0] ],
   #[ "V", "T", transpose pdl [  0,160, 0] ],

   #[ "I", "J", transpose pdl [  10, 10, 0] ],
   #[ "J", "I", transpose pdl [  10, 10, 0] ],
   #[ "I", "J", transpose pdl [-110, 10, 0] ],
   #[ "J", "I", transpose pdl [-110, 10, 0] ],

   #[ "M", "N", transpose pdl [  10, 10, 0] ],
   #[ "N", "M", transpose pdl [  10, 10, 0] ],
   #[ "M", "N", transpose pdl [ 110, 10, 0] ],
   #[ "N", "M", transpose pdl [ 110, 10, 0] ],

    [ "U", "I", transpose pdl [ -30, 20, 0] ],
    [ "V", "I", transpose pdl [ -30, 20, 0] ],
    [ "T", "I", transpose pdl [ -30, 20, 0] ],
); 

# &apply_sets;

# U(0.407769018436312)
# T(0.59223098156373)

# U(0.592230981563688)
# V(0.59223098156373)

# T(0.407769018436311)
# V(0.407769018436282)

# 0.408 0.592
my ($u_041, $u_059) = (0.407_769_018_436_312,
                       0.592_230_981_563_73);

my @a1a2a3 = qw/ q0  q1  n2 /;
my @array_c4c5c6 = (
    [qw/ p3 m3 m2 /],
    [qw/ p3 m3 m1 /],
    [qw/ m3 p3 m1 /],    # t3
    [qw/ m3 p3 m2 /],

    [qw/ p1 m1 m3 /],
    [qw/ p1 m1 m2 /],
    [qw/ m1 p1 m2 /],    # t7
    [qw/ m1 p1 m3 /],

    [qw/ p2 m2 m1 /],
    [qw/ p2 m2 m3 /],
    [qw/ m2 p2 m3 /],    # t11
    [qw/ m2 p2 m1 /],
);

&produce_transformation_matrices(\@a1a2a3, \@array_c4c5c6);
# &print_transformation_matrices;


### test routine(s)

### set grammar

my %elements;
my %tA_by_elements;
my %atom;

$elements{C} = new Set::Scalar qw/ U1 U2 U3 T4 T5 T6 V7 V8 V9 /;
$elements{A} = new Set::Scalar qw/ I J /;
$elements{B} = new Set::Scalar qw/ M N /;
$tA_by_elements{t3} = new Set::Scalar qw/ T5 T6 V8 V9 /;
$tA_by_elements{t7} = new Set::Scalar qw/ U1 U2 V7 V8 /;
$tA_by_elements{t11} = new Set::Scalar qw/ U2 U3 T4 T5 /;

print "C = $elements{C}\n";
print "tA{t3}  = $tA_by_elements{t3}\n";
print "tA{t7}  = $tA_by_elements{t7}\n";
print "tA{t11} = $tA_by_elements{t11}\n\n";

my $atom_index = "a";
$atom{$atom_index} = $elements{C};

while (my $this_transformation = each %tA_by_elements) {
    foreach my $this_atom (keys %atom) {
        my ($in, $out) = ($atom{$this_atom} * $tA_by_elements{$this_transformation},
                          $atom{$this_atom} - $tA_by_elements{$this_transformation});
        ($atom{$this_atom}, $atom{++$atom_index}) = ($in, $out) if $in->size and $out->size;
        $atom{$this_atom} = $in if $out->is_empty;
        $atom{$this_atom} = $out if $in->is_empty;
    }
}
undef $atom_index;

#while ( my $this = each %atom) {
#    print "atom $this = $atom{$this}\n";
#}

my %tA_by_atoms;
while (my $this_tA = each %tA_by_elements) {
    $tA_by_atoms{$this_tA} = Set::Scalar->new;
    while (my $this_atom = each %atom) {
        my $has_these = 1;
        $has_these &&= $tA_by_elements{$this_tA}->has($_) foreach @{ $atom{$this_atom} };
        $tA_by_atoms{$this_tA}->insert($this_atom) if $has_these;
    }
}


while ( my $this = each %tA_by_elements) {
    print "tA($this) = $tA_by_elements{$this}\n";
}
print "\n";

while ( my $this = each %tA_by_atoms) {
    print "tA($this) = $tA_by_atoms{$this}\n";
    my @atoms = @{ $tA_by_atoms{$this} };
    print "$_ has $atom{$_}\n" foreach @atoms;
}

my %maximal_elements;
my %points_on_element;
my %coords_on_element;

$maximal_elements {U} = [qw/ U1 U2 U3 /];
$points_on_element{U} = [qw/ p1 m2 m1 p2 /];
$coords_on_element{U} = [qw/ 0 $u_041 $u_059 1 /];

$maximal_elements {T} = [qw/ T4 T5 T6 /];
$points_on_element{T} = [qw/ p2 m3 m2 p3 /];
$coords_on_element{T} = [qw/ 0 $u_041 $u_059 1 /];

$maximal_elements {V} = [qw/ V7 V8 V9 /];
$points_on_element{V} = [qw/ p1 m3 m1 p3 /];
$coords_on_element{V} = [qw/ 0 $u_041 $u_059 1 /];

my $shape_C;
my %shape_tA;
my %shape_tB;



### main programme print $a, $c;

&gui;

my $canvas;

sub gui {
    my $mw = new MainWindow;
    $mw->title("Shape Grammars Synthesiser 2");
    $mw->geometry("+20+20");
    $mw->bind('<Double-Button-3>' => \&exit);
    $mw->bind('<Q>' => \&exit);

    &pull_down_menu($mw);
    my $t = &right_menu($mw);
    $canvas = &main_canvas($mw);
    &draw_things($canvas);
    #&draw_tB($canvas);

    MainLoop;
}

sub right_menu {
    my $mw = $_[0];
    my $right = $mw->Frame()->pack(
        -side => 'right',
        -anchor => 'n',
        -expand => '1',
        -fill => 'x',
    );

    $right->Label(
        -text => "Select a transformation",
    )->pack();

    my $list = $right->Scrolled('Listbox',
        -scrollbars => 'e',
        -selectmode => 'single',
    )->pack();

    # $list->insert('end', sort keys %transformation);
    $list->insert('end', qw/ t2 t4 t6 t8 t10 t12 /);

    my %on_off;

    $list->bind('<Double-Button-1>', [\&select_a_transformation, $list]);    

    my $mw_frame_buttons = $right->Frame()->pack(
        -side => 'right',
        -anchor => 'n',
        -expand => '1',
        -fill => 'x',
    );

    &toggle_buttons($mw_frame_buttons, \%on_off);
}

sub toggle_buttons {
    my $f = $_[0];
    my %on_off = %{$_[1]};
    my $b1 = $f->Checkbutton(
        -text => 'U',
        -variable => \$on_off{U},
    )->pack(-side => 'left');
}

sub select_a_transformation {
   my $item = $_[0];
   my $t = $item->get( $item->curselection() );
   #print "transformation $t is:\n";
   #print $transformation{$t};
   &draw_tB($canvas, $t);
}

sub draw_things {
    my $c = $_[0];
    &draw_C($c, "U", 0, 1, "U", 1, 'black');
    &draw_C($c, "T", 0, 1, "V", 1, 'black');
    &draw_C($c, "V", 0, 1, "T", 1, 'black');

    &draw_C($c, "U",      0, $u_041, "U1", 5, 'red');
    &draw_C($c, "U", $u_041, $u_059, "U2", 5, 'orange');
    &draw_C($c, "U", $u_059,      1, "U3", 5, 'green');

    &draw_C($c, "T",      0, $u_041, "T4", 5, 'green');
    &draw_C($c, "T", $u_041, $u_059, "T5", 5, 'pink');
    &draw_C($c, "T", $u_059,      1, "T6", 5, 'blue');

    &draw_C($c, "V",      0, $u_041, "V7", 5, 'red');
    &draw_C($c, "V", $u_041, $u_059, "V8", 5, 'yellow');
    &draw_C($c, "V", $u_059,      1, "V9", 5, 'blue');

    #&draw_C($c, "I", 0, 1, "I", 1, 'black');
    #&draw_C($c, "J", 0, 1, "J", 1, 'black');
    #&draw_C($c, "M", 0, 1, "M", 1, 'black');
    #&draw_C($c, "N", 0, 1, "N", 1, 'black');
}

sub draw_tB {
    my $c = $_[0];
    my $t = $_[1];

    #print "fdfdsfdfd transformation $t is:\n";
    #print $transformation{$t};

    my $MMM = &calc_tB_M("M", $t);
    my $NNN = &calc_tB_N("N", $t);
    
    &draw_C($c, $MMM, 0, 1, "M", 1, 'black');
    &draw_C($c, $NNN, 0, 1, "N", 1, 'black');
}

sub calc_tB_M {
    my ($C, $t) = @_;
    #print "curves === $C\n";
    #print "$t\n";
    my ($p0, $p1, $p2) = @{$arc{$C}};
    my $P0 = transpose pdl [ $coords{$p0}->at(0,0),
                             $coords{$p0}->at(0,1),
                             $coords{$p0}->at(0,2),
                             1];
    my $P1 = transpose pdl [ $coords{$p1}->at(0,0),
                             $coords{$p1}->at(0,1),
                             $coords{$p1}->at(0,2),
                             1];
    my $P2 = transpose pdl [ $coords{$p2}->at(0,0),
                             $coords{$p2}->at(0,1),
                             $coords{$p2}->at(0,2),
                             1];
    
    $P0 = $transformation{$t} x $P0;
    $P1 = $transformation{$t} x $P1;
    $P2 = $transformation{$t} x $P2;
    #print $P0, $P1, $P2;

    $coords{b0} = transpose pdl [ $P0->at(0,0),
                                  $P0->at(0,1),
                                  $P0->at(0,2) ];
    $coords{b1} = transpose pdl [ $P1->at(0,0),
                                  $P1->at(0,1),
                                  $P1->at(0,2) ];
    $coords{b2} = transpose pdl [ $P2->at(0,0),
                                  $P2->at(0,1),
                                  $P2->at(0,2) ];
    $arc{MM} = [qw/ b0 b1 b2 /];

    $weights{MM}   = [ &return_weights( @{$arc{MM}} ) ];
    $curvature{MM} = &return_curvature( @{$arc{MM}} );
    #print "arc MM : w1=$weights{MM}[1], 1/r=$curvature{MM}\n";

    return "MM";
}

sub calc_tB_N {
    my ($C, $t) = @_;
    #print "curves === $C\n";
    #print "$t\n";
    my ($p0, $p1, $p2) = @{$arc{$C}};
    my $P0 = transpose pdl [ $coords{$p0}->at(0,0),
                             $coords{$p0}->at(0,1),
                             $coords{$p0}->at(0,2),
                             1];
    my $P1 = transpose pdl [ $coords{$p1}->at(0,0),
                             $coords{$p1}->at(0,1),
                             $coords{$p1}->at(0,2),
                             1];
    my $P2 = transpose pdl [ $coords{$p2}->at(0,0),
                             $coords{$p2}->at(0,1),
                             $coords{$p2}->at(0,2),
                             1];
    
    $P0 = $transformation{$t} x $P0;
    $P1 = $transformation{$t} x $P1;
    $P2 = $transformation{$t} x $P2;
    #print $P0, $P1, $P2;

    $coords{b4} = transpose pdl [ $P0->at(0,0),
                                  $P0->at(0,1),
                                  $P0->at(0,2) ];
    $coords{b5} = transpose pdl [ $P1->at(0,0),
                                  $P1->at(0,1),
                                  $P1->at(0,2) ];
    $coords{b6} = transpose pdl [ $P2->at(0,0),
                                  $P2->at(0,1),
                                  $P2->at(0,2) ];
    $arc{NN} = [qw/ b4 b5 b6 /];

    $weights{NN}   = [ &return_weights( @{$arc{NN}} ) ];
    $curvature{NN} = &return_curvature( @{$arc{NN}} );
    #print "arc NN : w1=$weights{NN}[1], 1/r=$curvature{NN}\n";

    return "NN";
}

sub draw_C {
    my ($canvas, $C, $a, $b, $tag, $thickness, $colour) = @_;
    # print "$C\n";
    my @pairs_of_xy = ();
    my $u = $a;
    my $du =  ($b-$a)/20;
    for (my $i=0; $i<=20; $i++) {
        my $point = &C($C, $u);
        $u = $u + $du;
        $point *= 2;
        $point += transpose pdl [400, 300, 0];
        my ($x, $y) = ( $point->at(0,0), -$point->at(0,1) +300+500);
        push @pairs_of_xy, ($x, $y);
    };
    $canvas->createLine(@pairs_of_xy,
        -tag => $tag,
        -fill => $colour,
        -width => $thickness,
    );
}

sub pull_down_menu {
    my $mw = $_[0];
    my $pull_down = $mw->Frame()->pack(
        -side => 'top',
        -anchor => 'n',
        -expand => '1',
        -fill => 'x',
    );

    my $pd_menu_9 = $pull_down->Menubutton(
        -text => "Help",
        -tearoff => 0,
        -menuitems => [
            ['command' => "About",
             -command => \&mw_help_about,
            ],
        ],
    )->pack(-side => 'right');

    my $pd_menu_1 = $pull_down->Menubutton(
        -text => "File",
        -menuitems => [
            ['command' => "triangles",
             -command => \&mw_file_triangles,
            ],
            ['command' => "trefoil",
             -command => \&mw_file_trefoil,
            ],
            '-',
            ['command' => "Quit",
            -command => \&exit,
            -accelerator => "Ctrl-Q",
            ]
        ],
    )->pack(-side => 'left');    
}

sub main_canvas {
    my $mw = $_[0];
    my $canvas = $mw->Canvas(
        -width => 800,
        -height => 600,
        -cursor => 'tcross',
        -background => 'white',
    )->pack();
    return $canvas;
}

sub mw_help_about {
    print "about is here.\n";
}

sub mw_file_triangles {
    print "Open triangles\n";
}

sub mw_file_trefoil {
    print "Open trefoil\n";
}

exit 1;


### subroutines


sub print_transformation_matrices {
    while (my ($id, $matrix) = each %transformation) {
        print "transformation{$id} = $matrix";
    }
}

sub produce_transformation_matrices {
    my ($ref_a, $ref_array_c) = @_;
    my @a1a2a3       = @{$ref_a};
    my @array_c4c5c6 = @{$ref_array_c};
    # print "@a1a2a3\n";
    # print "@array_c4c5c6\n";

    my $i = 0;
    foreach my $item (@array_c4c5c6) {
        $i++;
        my @c4c5c6 = @{$item};
        my $name = "t" . $i;
        my $name_mirror = $name . "m";
        ($transformation{$name},
         $transformation{$name_mirror} )
            = &T_combined(@a1a2a3, @c4c5c6);
        # print "tranformation $name is defined by @c4c5c6\n\n";
        # print $name, "\n", $name_mirror, "\n";
        # print $transformation{$name}, $transformation{$name_mirror};
    };
    # the above produces (t1, t1m, ..., t12, t12m) -- 24 transformations
}

sub T_combined {
    my $a1 = $coords{$_[0]};
    my $a2 = $coords{$_[1]};
    my $a3 = $coords{$_[2]};
    my $c4 = $coords{$_[3]};
    my $c5 = $coords{$_[4]};
    my $c6 = $coords{$_[5]};

    my $T_a1_inv  = inv &T_translate($a1);
    my $T_c4      =     &T_translate($c4);

    my $A_uvw_inv = inv &T_rotate($a1, $a2, $a3);
    my $A_abc     =     &T_rotate($c4, $c5, $c6);
    
    my $S_s = &T_scale($a1, $a2, $c4, $c5);
    my $T_mirror = &T_mirror_zx;

    return ($T_c4 x $A_abc x $S_s             x $A_uvw_inv x $T_a1_inv,
            $T_c4 x $A_abc x $S_s x $T_mirror x $A_uvw_inv x $T_a1_inv);
}

sub T_mirror_zx {
    return pdl [ [1, 0, 0, 0],
                 [0,-1, 0, 0],
                 [0, 0, 1, 0],
                 [0, 0, 0, 1] ];
}

sub T_scale{
    my ($a1, $a2, $c4, $c5) = @_;
    my $s = &magnitude($c5 - $c4) / &magnitude($a2 - $a1);
    return pdl [ [$s,  0,  0, 0],
                 [ 0, $s,  0, 0],
                 [ 0,  0, $s, 0],
                 [ 0,  0,  0, 1] ];
}


sub T_translate {
    my $a1 = $_[0];
    return pdl [ [1, 0, 0, $a1->at(0,0)],
                 [0, 1, 0, $a1->at(0,1)],
                 [0, 0, 1, $a1->at(0,2)],
                 [0, 0, 0, 1           ] ];
}

sub T_rotate {
    my ($a1, $a2, $a3) = @_;
    my $a2a1 = $a2 - $a1;
    my $a3a1 = $a3 - $a1;
    my $unit_vector_u = norm squeeze $a2a1;
    my $unit_vector_w = norm crossp $unit_vector_u, squeeze $a3a1;
    my $unit_vector_v = norm crossp $unit_vector_w, $unit_vector_u;
    return pdl [$unit_vector_u,
                $unit_vector_v,
                $unit_vector_w,
                [0, 0, 0, 1] ];
}

sub apply_sets {
    foreach (@sets) {
        my ($U, $T, $P) = @{$_};
        my ($point, $u, $t) = &find_intersection($P, $U, $T);
        print "Intersection point found ...\n";
        print "\t$U($u)\n";
        print "\t$T($t)";
        print &C($U, $u);
        print "--- * --- * ---\n";
    }
}

sub find_intersection {
    my ($P, $U, $T) = @_;
    my $u = my $initial_u = &closest_u($P, $U);
    my $t = my $initial_t = &closest_u($P, $T);
    my ($P_u, $P_t);
    $P_t = $P;
    my $d = 'Inf';
    my $last_d;
    my $i = 0;
    LOOP: while ($d > $epsilon1 * 0.01) {
        $i++;
        if ($i > 100) {
            print "=== * === * ===\n";
            print "vvv No intersection found vvv\n";
            print "initial $U($initial_u)\n";
            print "current $U($u)= $P_u";
            print "--- * --- * ---\n";
            print "initial $T($initial_t)\n";
            print "current $T($t)= $P_t";
            print "--- * --- * ---\n";
            print "distance=$d\n";
            print "delta_distance=", $last_d-$d, "\n";
            print "=== * === * ===\n";
            print "^^^ looped too much time in subroutine ", (caller(0))[3], " ^^^\n";
            last LOOP;
        }

        $u = &point_projection($P_t, $U, $u);
        $P_u = &C($U, $u);

        $t = &point_projection($P_u, $T, $t);
        $P_t = &C($T, $t);

        $last_d = $d;
        $d = &distance($P_u, $P_t);
    }

    print "loop number $i\n";
    print "distance = $d\n";
    # print "$U($u)\n";
    # print "$T($t)\n";
    # print "[ ", $P_t->at(0,0), " ]\n";
    # print "[ ", $P_t->at(0,1), " ]\n";
    # print "[ ", $P_t->at(0,2), " ]\n";
    return ($P_t, $u, $t);
}

sub closest_u {
    my ($P, $U) = @_;
    my $u=0;
    my $d='Inf';
    for (my $i=0; $i<=10; $i++) {
        my $P_u = &C($U, 0.1 * $i);
        my $new_d = &distance($P, $P_u);
        if ($new_d < $d) {
            $d = $new_d;
            $u = 0.1 * $i;
        }
    }
    return $u;
}

sub point_projection {
    my ($P, $arc, $u) = @_;
    my $Q = &C($arc, $u);
    my $d = &distance($P, $Q);
    while ($d > $epsilon1 * 0.01) {
        my $C1 = &C_derv($arc, $u);
        my $QP = $Q - $P;
        # dot product in PDL = sum($A*$B);
        my $f  = sum( $C1 * $QP );
        my $f1 = sum( &C_derv2($arc, $u) * $QP ) +
                 sum( $C1 ** 2 );
        my $u_last = $u;
        $u -= $f / $f1;
        $Q = &C($arc, $u);
        $d = &distance($P, $Q);
        my $C1_magnitude = &magnitude($C1);
        last if abs ($u-$u_last) * $C1_magnitude < $epsilon1 * 0.001 or
                $d < $epsilon1 * 0.001 or
                abs($f) / ($C1_magnitude*$d) < $epsilon2;
    }
    $u=0 if $u < 0;
    $u=1 if $u > 1;
    return $u;
}

sub magnitude {
    my $d = $_[0];
    my ($x, $y, $z) = ( $d->at(0,0), $d->at(0,1), $d->at(0,2) );
    return sqrt $x**2 + $y**2 + $z**2;
}

sub distance {
    my ($P, $Q) = @_;
    my $PQ = $P - $Q;
    my ($dx, $dy, $dz) = ( $PQ->at(0,0), $PQ->at(0,1), $PQ->at(0,2) );
    return sqrt $dx**2 + $dy**2 + $dz**2;
}

sub C {
    my ($arc, $u) = @_;
    my $A = &A($arc, $u);
    my $w = &w($arc, $u);
    return $A / $w;
}

sub C_derv {
    my ($arc, $u) = @_;
    my $A  = &A     ($arc, $u);
    my $A1 = &A_derv($arc, $u);
    my $w  = &w     ($arc, $u);
    my $w1 = &w_derv($arc, $u);
    return ($w*$A1 - $w1*$A) / $w**2;
}

sub C_derv2 {
    my ($arc, $u) = @_;
    my $A  = &A      ($arc, $u);
    my $A1 = &A_derv ($arc, $u);
    my $A2 = &A_derv2($arc, $u);
    my $w  = &w      ($arc, $u);
    my $w1 = &w_derv ($arc, $u);
    my $w2 = &w_derv2($arc, $u);
    return (    $w**2          * $A2
            -   $w    * $w2    * $A
            - 2*$w    * $w1    * $A1
            + 2       * $w1**2 * $A) / $w**3;
}

sub A {
    my ($arc, $u) = @_;
    my ($p0, $p1, $p2) = @{ $arc{$arc} };
    my $P0 = $coords{$p0};
    my $P1 = $coords{$p1};
    my $P2 = $coords{$p2};
    my ($w0, $w1, $w2) = @{ $weights{$arc}};
    return (1-$u)**2 * $P0 + 2*$u*(1-$u) * $w1 * $P1 + $u**2 * $P2;
}

sub A_derv {
    my ($arc, $u) = @_;
    my ($p0, $p1, $p2) = @{ $arc{$arc} };
    my $P0 = $coords{$p0};
    my $P1 = $coords{$p1};
    my $P2 = $coords{$p2};
    my ($w0, $w1, $w2) = @{ $weights{$arc}};
    return (-2 + 2*$u) * $P0 + (2 - 4*$u) * $w1 * $P1 + 2*$u * $P2;
}

sub A_derv2 {
    my ($arc, $u) = @_;
    my ($p0, $p1, $p2) = @{ $arc{$arc} };
    my $P0 = $coords{$p0};
    my $P1 = $coords{$p1};
    my $P2 = $coords{$p2};
    my ($w0, $w1, $w2) = @{ $weights{$arc}};
    return 2 * $P0 - 4*$w1 * $P1 + 2 * $P2;
}

sub w {
    my ($arc, $u) = @_;
    my ($w0, $w1, $w2) = @{ $weights{$arc}};
    return (1-$u)**2 + 2*$u*(1-$u) * $w1 + $u**2;
}

sub w_derv {
    my ($arc, $u) = @_;
    my ($w0, $w1, $w2) = @{ $weights{$arc}};
    return -2 * (1-$w1) + 4 * (1-$w1) * $u;
}

sub w_derv2 {
    my ($arc, $u) = @_;
    my ($w0, $w1, $w2) = @{ $weights{$arc}};
    return 4 * (1-$w1);
}

sub return_weights {
    my ($B, $A, $C) = @_;
    my $a = &distance( $coords{$B}, $coords{$C} );
    my $b = &distance( $coords{$C}, $coords{$A} );
    my $c = &distance( $coords{$A}, $coords{$B} );
    my $angle_A = acos( ($b**2 + $c**2 - $a**2) / (2 * $b * $c) );
    my $angle_theta = $pi/2 - $angle_A/2;
    # my $kappa = 1 / ($c * tan ($angle_A/2));
    return (1, cos $angle_theta, 1);
}

sub return_curvature {
    my ($B, $A, $C) = @_;
    my $a = &distance( $coords{$B}, $coords{$C} );
    my $b = &distance( $coords{$C}, $coords{$A} );
    my $c = &distance( $coords{$A}, $coords{$B} );
    my $angle_A = acos( ($b**2 + $c**2 - $a**2) / (2 * $b * $c) );
    my $angle_theta = $pi/2 - $angle_A/2;
    my $kappa = 1 / ($c * tan ($angle_A/2));
    return $kappa;
}

sub trefoil_control_polygon {
    my $d1 = sqrt( 80**2 - 55**2 );
    my $d2 = 55 / tan( $pi/3 );
    my $h1 = sqrt( $d2**2 + 55**2 );
    my $A = asin( sin($pi/3) * $h1 / 80 );
    my $B = $pi - $A - $pi/3;
    my $h2 = 80 * sin($B) / sin($pi/3);
    my $x1 = $h2 * sin($pi/3);
    # print "x1 = $x1\n";
    my $y1 = sqrt( 80**2 - $x1**2); 
    # print "y1 = $y1\n";
    my $d3 = $h2 * cos($pi/3);
    my $y2 = $d1 + $d2 + $d3 + $y1;
    # print "y2 = $y2\n";
    my $c3 = 80**2 / $y1; 
    # print "c3 = $c3\n";
    my $x3 = ($c3-$d3-$y1) * cos($pi/6);
    my $y3 =-($c3-$d3-$y1) * sin($pi/6) + $d3 + $y1;
    # print "x3 = $x3\n";
    # print "y3 = $y3\n";
    # print "0, 0\n";                       # centre of arc a3
    # print "55, ", $d2 + $d3 + $y1, "\n";  # centre of arc a2
    return ( $x1, $y1, $y2, $c3, $x3, $y3, 55, $d2 + $d3 + $y1);
}

sub trefoil_rule_polygon {
    my $x = $_[0] / 2;
    my $r = 80;
    my $adj = sqrt $r**2 - $x**2;
    my $H = $r**2 / $adj;
    my $y = $H - $adj; 

    # print "$_[0], $x, $y\n";

    $coords{q0} = transpose pdl [     0,  0, 0];
    $coords{q1} = transpose pdl [-$_[0],  0, 0];
    $coords{q2} = transpose pdl [ $_[0],  0, 0];

    $coords{n1} = transpose pdl [  -$x,  $y, 0];
    $coords{n2} = transpose pdl [  -$x, -$y, 0];
    $coords{n3} = transpose pdl [   $x,  $y, 0];
    $coords{n4} = transpose pdl [   $x, -$y, 0];

    $arc{I} = [qw/ q0 n1 q1 /];  # lhs of a rule
    $arc{J} = [qw/ q0 n2 q1 /];
    $arc{M} = [qw/ q0 n3 q2 /];  # rhs of a rule
    $arc{N} = [qw/ q0 n4 q2 /];
}

__END__