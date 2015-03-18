#!/usr/bin/perl
my $Nbins = 512;
my $inputbase = "/data3/sinham/ZA/512/z300/high_int_accuracy/run3/";
my $snapshotbase = "za_512_z300_PM_RUN3";
my $outputbase = "/data3/sinham/ZA/512/z300/high_int_accuracy/run3/";

for ($count = 0; $count <= 61; $count++) {
	system("./makehist_density",$Nbins,$inputbase,$snapshotbase,$count,$outputbase);
	exit 1 if ( $? != 0 );
}

