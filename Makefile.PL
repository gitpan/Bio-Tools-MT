use ExtUtils::MakeMaker;


WriteMakefile(
    'NAME'		=> 'Bio::Tools::MT',
    'VERSION_FROM'	=> 'MT.pm', 
    'PREREQ_PM'		=> {}, 
    ($] >= 5.005 ?    
      (ABSTRACT_FROM => 'MT.pm', 
       AUTHOR     => 'xern <xern@cpan.org>') : ()),
    'LIBS'		=> [''], 
    'DEFINE'		=> '', 
	
    'INC'		=> '', 
	
    
);
