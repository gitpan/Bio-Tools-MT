package Bio::Tools::MT;

use strict;
our $VERSION = '0.01';

use Exporter;
our @ISA = qw/Exporter/;
our @EXPORT_OK = qw/weight/;

sub new {
    my $pkg = shift;
    my $arg = { @_ };
    my $seq = join q//, (uc($arg->{sequence}) =~ /[ATCG]/go);

    bless {

	sequence        => ( $seq || die "No sequence is given\n"),

	kna_conc        => ($arg->{kna_conc} || 50),
	mg_conc         => ($arg->{mg_conc}  || 0),
	salt_conc       => ($arg->{kna_conc} || 50),

	primer_conc     => ($arg->{primer_conc} || 200),
	target_conc     => ($arg->{target_conc} || 0),

	has_pg          => ($arg->{has_pg} || 0),

    }, $pkg;
}


sub basecount {
    my $seqref = shift;
    my %cnt;
    $cnt{$&}++ while $$seqref =~ /[ATCG]/go;
    \%cnt;
}

# is self-complementary?
sub is_sc {
    my $seqref = shift;
    my %c = qw/A T T A C G G C/;
    my $newseq = $$seqref;
    $newseq =~ s/./$c{$&}/go;
    1 if $newseq eq $$seqref;
}


# ----------------------------------------------------------------------
sub calc_mt_0 {
    my $seqref = shift;
    my $cnt = shift;

    my $len = length $$seqref;

    return
	sprintf "%.0f",
	($len<14) ?
	4 * ($cnt->{G} + $cnt->{C}) + 2 * ($cnt->{A} + $cnt->{T}) :
	     64.9 + 41 * ($cnt->{G} + $cnt->{C} - 16.4) / $len;
}

# ----------------------------------------------------------------------
sub calc_mt_1 {
    my $seqref = shift;
    my $cnt = shift;
    my $kna_conc = shift;

    my $len = length $$seqref;

    return
	sprintf "%.0f",
	81.5 + 16.6 * (log($kna_conc*1E-3)/log(10)) + 0.41*(100*($cnt->{G} + $cnt->{C})/$len) - 675/$len;
}


# ----------------------------------------------------------------------

my %dhs = (
	   GG => [ qw/8    19.9/ ],
	   GA => [ qw/8.2  22.2/ ],
	   GT => [ qw/8.4  22.4/ ],
	   GC => [ qw/10.6 27.2/ ],
	   AG => [ qw/7.8  21/   ],
	   AA => [ qw/7.9  22.2/ ],
	   AT => [ qw/7.2  20.4/ ],
	   AC => [ qw/8.4  22.4/ ],
	   TG => [ qw/8.5  22.7/ ],
	   TA => [ qw/7.2  21.3/ ],
	   TT => [ qw/7.9  22.2/ ],
	   TC => [ qw/8.2  22.2/ ],
	   CG => [ qw/10.6 27.2/ ],
	   CA => [ qw/8.5  22.7/ ],
	   CT => [ qw/7.8  21/   ],
	   CC => [ qw/8    19.9/ ],
	   );


# helix initiation corrections from Santalucia 1998 & Allawi & Santalucia 1997
sub terminalcorrections {
    my $seqref = shift;

    my $deltah=0;
    my $deltas=0;
    my $first = substr($$seqref, 0, 1);
    my $last = substr($$seqref, length($$seqref)-1, 1);

    if( $first =~ /[GC]/o ){
	$deltah+=0.1;
	$deltas+=-2.8;
    }
    elsif( $first =~ /[AT]/o ){
	$deltah+=2.3;
	$deltas+=4.1;
    }
    
    if ( $last =~ /[GC]/o ){
	$deltah+=0.1;
	$deltas+=-2.8;
    }
    elsif( $last =~ /[AT]/o ){
	$deltah+=2.3;
	$deltas+=4.1;
    }
    return ($deltah, $deltas);
}

# changes to ds dependant on the salt concentration & sequence length
sub saltcorrections {
    my ($seqref, $param) = @_;

    # convert to moles and then adjust for greater stabilizing effects
    # of Mg compared to Na or K. See von Ahsen et al 1999

    my $deltas;
    $param->{salt_conc} = $param->{kna_conc} * 1E-3;
    $param->{salt_conc} += ($param->{mg_conc} * 1E-3 * 140);

    # This comes from von Ahsen et al 1999

    # Delta S
    return 0.368 * (length($$seqref)-1)* log($param->{salt_conc});#/log(10);
}


# base stacking calculations. 
sub calc_mt_2 {
    my ($seqref, $param) = @_;
    my ($ds, $dh) = (0, 0);

    $ds += saltcorrections($seqref, $param);

    my @t = terminalcorrections($seqref);
    $dh+=$t[0];
    $ds+=$t[1];

    for(my $i = 0 ; $i< length($$seqref)-1 ; $i++){
	$dh-=$dhs{substr($$seqref, $i, 2)}->[0];
	$ds-=$dhs{substr($$seqref, $i, 2)}->[1];
    }


    if( is_sc($seqref) ){
	sprintf "%.0f", (-273.15+ 1000 * $dh / ($ds+1.987*log(($param->{primer_conc}/4)*1E-9)))
    }
    else {
	sprintf "%.0f", ( -273.15 + (1000*$dh)/($ds+1.987*(log(($param->{primer_conc} - $param->{target_conc})*1E-9 / 2))));
    }

}

# ----------------------------------------------------------------------

sub result {
    my $pkg = shift;
    my $ret;

    # GC ratio
    my $cnt = basecount(\${$pkg}{sequence});

    $ret->{gc_ratio} = sprintf "%.2f", ($cnt->{G} + $cnt->{C})/length($pkg->{sequence});
    $ret->{length} = length($pkg->{sequence});

    $ret = {
	%$ret,
	%{$pkg->temperature($cnt)},
	weight => $pkg->weight(),
    };
}

sub temperature {
    my ($pkg, $cnt) = @_;
    my %ret;

    $cnt = basecount(\${$pkg}{sequence}, $cnt) unless ref $cnt;

    # model 0
    $ret{model_0} = calc_mt_0(\${$pkg}{sequence}, $cnt);

    # model 1
    $ret{model_1} = calc_mt_1(\${$pkg}{sequence}, $cnt, $pkg->{kna_conc});

    # model 2
    $ret{model_2} = calc_mt_2(\${$pkg}{sequence}, $pkg);

    \%ret;
}

sub weight {
    my $pkg;
    my $has_pg = 0;
    my $seq;

    if(ref($_[0]) eq 'Bio::Tools::MT'){
	$pkg = shift;
	$seq = $pkg->{sequence};
	$has_pg = $pkg->{has_pg};
    }
    else{
	$seq = shift;
	$has_pg = shift;
    }
    my $cnt = basecount(\$seq);


    # -61 => (-79 for missing phosphate (PO3) at the 5'end, +1 for the hydrogen
    # that replaces it & +17 for a 3'hydroxyl)  
    return sprintf "%.0f", -61 + (330.2 * $cnt->{G}) + (314.2 * $cnt->{A}) + (305.2 * $cnt->{T}) + (290.2 * $cnt->{C}) + ($has_pg ? +78 : 0);
    
}





1;



__END__
# Below is stub documentation for your module. You better edit it!

=head1 NAME

Bio::Tools::MT - Calculating melting temperature of oligos

=head1 SYNOPSIS

  use Bio::Tools::MT;
  my $mt = Bio::Tools::MT->new(
			       sequence => $ATCG,
			       );
  print $mt->temperature;
  print $mt->weight;


=head1 MELTING TEMPERATURE

There are various ways to calculate the melting temperature (Tm) of an oligo, and all of these will yield different results. This module tries to bundle available methods together for your use. 

The formula the module knows are as follows,

=head2 Simplest formula

Tm = 4 * (number of G's and C's in the primer) + 2 * (number of A's and T's in the primer)

but this only works for oligos < 14 bases and the reaction happens in presence of 50mM monovalent cations. For longer ones, another approximation is adopted.

Tm =  64.9 + 41 * (number of G's and C's in the primer - 16.4)/N


=head2 Salt-adjusted Calculation

Another formula takes into consideration the salt concentration of the reaction.

Tm =  81.5  +  16.6 * (log10([Na+] + [K+]))  +  0.41 * (%GC) - 675/N

=head2 Base-stacking Calculation

Another more advanced one takes also base stacking parameters.

If the primer is a non-self-complementary oligo and the primer concentration is greater than the target concentration, then

Tm = ( delta(H) / ( delta(S) + R * ln([primer]/2) ) ) - 273.15

I<delta(H)> is the enthalpy of base stacking interactions adjusted for helix initiation factors.

I<delta(S)> is the entropy of base stacking adjusted for helix initiation factors (6,7) and for the contributions of salts(a)  to the entropy of the system (6).

I<R> is the universal gas constant.


For self-complementary oligos, the denominator becomes delta(S) + R*ln([primer]/4). And if the primer concentration and the target concentration are almost equal, the denominator becomes delta(S) + R ln([primer] - [target]/2).

For polymers, the situation gets more complex and the salt effects on polymers are quite different. See B<SantaLucia, J. (1998) Proc. Nat. Acad. Sci. USA 95, 1460.>

=head1 MOLECULAR WEIGHT

In addition to melting temperature calculation, this module also calculates the molecular weight of a given sequence.

Molecular weight = (330.2 * number of G's) + (314.2 * number of A's) + (305.2 * number of T's) + (290.2 * number of C's), but the weight must be adjusted by +78 if a 5'phosphate is present. Default is unadjusted.

=head1 INTERFACE

=head2 new

Takes the background settings.

    my $mt = Bio::Tools::MT->new(
				 sequence    => $ATCG,
				 kna_conc    => 50,  # 50 mM
				 primer_conc => 200, # 200 nanomolar
                                 target_cont => 190, # 190 nanomolar
				 mg_conc     => 1.5, # 1.5 mM
				 has_pg      => 1,   # has a phosphate group
				 );

Characters of sequence other than ATCG will be ignored.

I<kna_conc> is for the combined concentration of K+ and Na+. Default is 50 (mH)

I<primer_conc> for the primer concentration in the reaction. (N.B. Primer concentration only affects the base-stacking calculations.) Default is 200 nonamolar.

I<target_conc> for the target concentration. It is set if it is close to the primer concentration, and when calculating temperature using the third model, the denominator changes a little.

I<mg_conc> for the Mg(2+) concentration. This only affects the base-stacking calculations.

I<has_pg> is set to 1 if the oligo has a 5'-phosphate group. (N.B. This affects only the molecular weight, not the temperature. The weight is adjusted by +78 if a 5'phosphate is present.) Default is 0.


=head2 result

    $mt->result();

Returns an anonymous hash containing the statistics about the sequence, such as temperature, molecular weight, gc-ratio, etc.

=head2 temperature

Outputs the predicted melting temperature.

    use Data::Dumper;
    print $mt->temperature();

In the returned anonymous hash, I<model> stands for the above 3 temperature calculation formulae: 0 for the simplest, 1 for the salt-adjusted, and 2 for the advanced one.


=head2 weight

Returns its molecular weight

    print $mt->weight();

Also, you have a non-OOish way.

    use Bio::Tools::MT qw/weight/;

    print weight($ATCG);

    print weight($ATCG, 1);  # the one follows if a 5'phosphate group is present


=head1 COPYRIGHT

xern E<lt>xern@cpan.orgE<gt>

This module is free software; you can redistribute it or modify it under the same terms as Perl itself.

=cut
