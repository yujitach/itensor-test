for($gam=0.3;$gam<0.38;$gam+=0.005){
	for($a2=-0.063;$a2<-0.062;$a2+=0.001){
		$filename=sprintf("%.10f %.10f .compute",$gam,$a2);
		system(qq(touch "data/$filename"));
	}
}