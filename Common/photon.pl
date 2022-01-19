use strict;
require JSON::PP;
use JSON::PP;


# X-RAY DATA BOOKLET
# Center for X-ray Optics and Advanced Light Source

my %Kedges = ('H', 13.6, 'He', 24.6, 'Li', 54.7, 'Be', 111.5, 'B', 188, 'C', 284.2, 
              'N', 409.9, 'O', 543.1, 'F', 696.7, 'Ne', 870.2, 'Na', 1070.8, 'Mg', 1303, 
              'Al', 1559, 'Si', 1839, 'P', 2145.5, 'S', 2472, 'Cl', 2822.4, 
              'Ar', 3205.9, 'K', 3608.4, 'Ca', 4038.5, 'Sc', 4492, 'Ti', 4966, 
              'V', 5465, 'Cr', 5989, 'Mn', 6539, 'Fe', 7112, 'Co', 7709, 
              'Ni', 8333, 'Cu', 8979, 'Zn', 9659, 'Ga', 10367, 'Ge', 11103,
              'As', 11867, 'Se', 12658, 'Br', 13474, 'Kr', 14326, 'Rb', 15200,
              'Sr', 16105, 'Y', 17038, 'Zr', 17998, 'Nb', 18986, 'Mo', 20000, 
              'Tc', 21044, 'Ru', 22117, 'Rh', 23220, 'Pd', 24350, 'Ag', 25514, 
              'Cd', 26711, 'In', 27940, 'Sn', 29200, 'Sb', 30491, 'Te', 31814, 
              'I', 33169, 'Xe', 34561, 'Cs', 35985, 'Ba', 37441, 'La', 38925, 
              'Ce', 40443, 'Pr', 41991, 'Nd', 43569, 'Pm', 45184, 'Sm', 46834, 
              'Eu', 48519, 'Gd', 50239, 'Tb', 51996, 'Dy', 53789, 'Ho', 55618, 
              'Er', 57486, 'Tm', 59390, 'Yb', 61332, 'Lu', 63314, 'Hf', 65351, 
              'Ta', 67416, 'W', 69525, 'Re', 71676, 'Os', 73871, 'Ir', 76111, 
              'Pt', 78395, 'Au', 80725, 'Hg', 83102, 'Tl', 85530, 'Pb', 88005, 
              'Bi', 90526, 'Po', 93105, 'At', 95730, 'Rn', 98404, 'Fr',101137, 
              'Ra',103922, 'Ac',106755, 'Th',109651, 'Pa',112601, 'U',115606 );

my %L2edges = ('F', 48.5, 'Ne',63.5, 'Na',88.6, 'Mg',17.8, 'Al',49.7, 'Si',89,
               'P', 136, 'S', 163.6, 'Cl',202, 'Ar',250.6, 'K', 297.3, 'Ca',349.7,
               'Sc',403.6, 'Ti',460.2, 'V', 519.8, 'Cr',583.8, 'Mn',649.9, 'Fe',719.9,
               'Co',793.2, 'Ni',870, 'Cu',952.3, 'Zn', 1044.9, 'Ga', 1143.2, 'Ge', 1248.1,
               'As', 1359.1, 'Se', 1474.3, 'Br', 1596, 'Kr', 1730.9, 'Rb', 1864,
               'Sr', 2007, 'Y',  2156, 'Zr', 2307, 'Nb', 2465, 'Mo', 2625,
               'Tc', 2793, 'Ru', 2967, 'Rh', 3146, 'Pd', 3330, 'Ag', 3524,
               'Cd', 3727, 'In', 3938, 'Sn', 4156, 'Sb', 4380, 'Te', 4612,
               'I',  4852, 'Xe', 5107, 'Cs', 5359, 'Ba', 5624, 'La', 5891,
               'Ce', 6164, 'Pr', 6440, 'Nd', 6722, 'Pm', 7013, 'Sm', 7312,
               'Eu', 7617, 'Gd', 7930, 'Tb', 8252, 'Dy', 8581, 'Ho', 8918,
               'Er', 9264, 'Tm', 9617, 'Yb', 9978, 'Lu', 10349, 'Hf', 10739,
               'Ta', 11136, 'W',  11544, 'Re', 11959, 'Os', 12385, 'Ir', 12824,
               'Pt', 13273, 'Au', 13734, 'Hg', 14209, 'Tl', 14698, 'Pb', 15200,
               'Bi', 15711, 'Po', 16244, 'At', 16785, 'Rn', 17337, 'Fr', 17907,
               'Ra', 18484, 'Ac', 19083, 'Th', 19693, 'Pa', 20314, 'U' , 20948);

my %L3edges = ('Ne',21.6, 'Na',30.5, 'Mg',49.2, 'Al',72.5, 'Si',99.2, 'P', 135, 'S', 162.5, 
              'Cl',200, 'Ar',248.4, 'K', 294.6, 'Ca',346.2, 'Sc',398.7, 'Ti',453.8,
              'V', 512.1, 'Cr',574.1, 'Mn',638.7, 'Fe',706.8, 'Co',778.1, 'Ni',852.7,
              'Cu',932.7, 'Zn', 1021.8, 'Ga', 1116.4, 'Ge', 1217, 'As', 1323.6,
              'Se', 1433.9, 'Br', 1550, 'Kr', 1678.4, 'Rb', 1804, 'Sr', 1940,
              'Y',  2080, 'Zr', 2223, 'Nb', 2371, 'Mo', 2520, 'Tc', 2677,
              'Ru', 2838, 'Rh', 3004, 'Pd', 3173, 'Ag', 3351, 'Cd', 3538,
              'In', 3730, 'Sn', 3929, 'Sb', 4132, 'Te', 4341, 'I',  4557,
              'Xe', 4786, 'Cs', 5012, 'Ba', 5247, 'La', 5483, 'Ce', 5723,
              'Pr', 5964, 'Nd', 6208, 'Pm', 6459, 'Sm', 6716, 'Eu', 6977,
              'Gd', 7243, 'Tb', 7510, 'Dy', 7790, 'Ho', 8071, 'Er', 8358,
              'Tm', 8648, 'Yb', 8944, 'Lu', 9244, 'Hf', 9561, 'Ta', 9881,
              'W',  10207, 'Re', 10535, 'Os', 10871, 'Ir', 11215, 'Pt', 11564,
              'Au', 11919, 'Hg', 12284, 'Tl', 12658, 'Pb', 13035, 'Bi', 13419,
              'Po', 13814, 'At', 14214, 'Rn', 14619, 'Fr', 15031, 'Ra', 15444,
              'Ac', 15871, 'Th', 16300, 'Pa', 16733, 'U' , 17166) ;



  my @z2symb =          ( '', 'H' , 'He', 'Li', 'Be', 'B' , 'C' , 'N' ,
      'O' , 'F' , 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P' , 'S' , 'Cl', 'Ar',
      'K' , 'Ca', 'Sc', 'Ti', 'V' , 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu',
      'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Rb', 'Sr', 'Y' , 'Zr',
      'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn', 'Sb',
      'Te', 'I' , 'Xe', 'Cs', 'Ba', 'La', 'Ce', 'Pr', 'Nd', 'Pm',
      'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu', 'Hf',
      'Ta', 'W' , 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb', 'Bi',
      'Po', 'At', 'Rn', 'Fr', 'Ra', 'Ac', 'Th', 'Pa', 'U' , 'Np', 'Pu',
      'Am', 'Cm', 'Bk', 'Cf', 'Es', 'Fm', 'Md', 'No', 'Lr', 'Rf', 'Db',
      'Sg', 'Bh', 'Hs', 'Mt' );


my $oceanData;
my $json = JSON::PP->new;
my $dataFile = "postDefaultsOceanDatafile";
if( open my $in, "<", $dataFile ) {
  local $/ = undef;
  $oceanData = $json->decode(<$in>);
  close($in);
} else {
  die "Failed to open $dataFile\n$!";
}


exit 0 if( scalar @{$oceanData->{'calc'}->{'edges'}}< 1 );


$oceanData->{'calc'}->{'edges'}[0] =~ m/(\d+)\s+(\d+)\s+(\d+)/ or die "Malformed calc->edges\n";
my $site = $1;
my $n = $2;
my $l = $3;

my $typat = $oceanData->{'structure'}->{'typat'}[$site-1];
my $znucl = $oceanData->{'structure'}->{'znucl'}[$typat-1];

my $zsymb = $z2symb[$znucl];
#print "$site $typat $znucl $zsymb\n";

my $energy = 0;
if( $n == 1 ) {
  $energy = $Kedges{ $zsymb };
} elsif( $n == 2 ) {
  if( $l == 1 ) {
    $energy = $L2edges{ $zsymb } + 2*$L3edges{ $zsymb };
    $energy/= 3;
  }
} else {
  print "WARNING Unsupported edge!\n";
}

#print "$energy\n";

my @dir = ( [ 1, 0, 0 ], [ 0, 1, 0 ], [0, 0, 1] );

my $quad = 0;
$quad = 1 if( $energy > 4000 ) ;

my $nphoton = 0;
for( my $i = 0; $i < 3; $i++ ) {
  for( my $j = 0; $j<3; $j++ ) {
    next if( $j == $i );
    $nphoton ++;
    my $fileName = sprintf "default_photon%i", $nphoton;
    open OUT, ">", $fileName or die;
    if( $quad == 1 ) {
      print OUT "quad\n";
    } else {
      print OUT "dipole\n";
    }
    printf OUT "cartesian %i %i %i\nend\n", $dir[$i][0], $dir[$i][1], $dir[$i][2];
    printf OUT "cartesian %i %i %i\nend\n", $dir[$j][0], $dir[$j][1], $dir[$j][2];
    printf OUT "%.1f\n", $energy;
    close OUT;
    last unless( $quad == 1 );
  }
}
