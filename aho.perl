#$gamstart=$ARGV[0]+0.0025;
#$astart=$ARGV[1]+0.00025;

gen(0.25,-0.05);
gen(0.30,-0.07);
gen(0.35,-0.09);
gen(0.40,-0.11);

sub gen{
	my ($gs,$as)=@_;
	generate($gs,$as);
#	generate($gs+0.0025,$as+0.00025);
}
sub generate{
	($gstart,$astart)=@_;
for($gam=$gstart;$gam<$gstart+0.2;$gam+=0.005){
	#$astart=-0.4*$gam+0.05;
	for($a2=$astart;$a2<$astart+0.02;$a2+=0.0005){
		$filename=sprintf("%.10f %.10f .compute",$gam,$a2);
		system(qq(touch "data/$filename"));
	}
}
}
