# c.f. http://www.promega.com/biomath/calc11.htm

use Test;
BEGIN { plan tests => 8 };


use Bio::Tools::MT qw/weight/;
ok(1);


# ----------------------------------------------------------------------
$seq = 'ACACAGTAGTGTCGTG';
my $mt = Bio::Tools::MT->new(
			     sequence    => $seq,
			     kna_conc    => 50,  # 50 mM
			     primer_conc => 200, # 200 nanomolar
			     );
ok($mt->result()->{'length'}, 16);
ok($mt->result()->{'model_2'}, 48);
ok($mt->temperature()->{'model_1'}, 38);
ok($mt->weight, 4938);

# ----------------------------------------------------------------------
$seq = 'AGCTTGCGTG';

my $mt = Bio::Tools::MT->new(
			     sequence    => $seq,
			     kna_conc    => 50,  # 50 mM
			     primer_conc => 200, # 200 nanomolar
			     mg_conc     => 1.5, # 1.5 mM
			     has_pg      => 1,   # has a phosphate group
			     );
ok($mt->result()->{'model_0'}, 32);
ok($mt->temperature()->{'model_2'}, 42);
ok($mt->weight, weight($seq, 1));

