for($gam=0.25;$gam<0.35;$gam+=0.0025){
#$a2=-0.02;
	for($a2=-0.05;$a2<-0.03;$a2+=0.00025){
		$filename=sprintf("%.10f %.10f .compute",$gam,$a2);
		system(qq(touch "data/$filename"));
	}
}
