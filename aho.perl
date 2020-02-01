for($gam=0.2;$gam<0.4;$gam+=0.005){
#$a2=-0.02;
	for($a2=-0.07;$a2<-0.06;$a2+=0.0005){
		$filename=sprintf("%.10f %.10f .compute",$gam,$a2);
		system(qq(touch "data/$filename"));
	}
}
