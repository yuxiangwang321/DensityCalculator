# DensityCalculator
# file_1 file_2: need these two files for calculation
# type : mass_xy, mass_yz, mass_xz, mass_x, mass_y, mass_z are for mass density caculation;
#        number_xy, number_yz, number_xz, number_x, number_y, number_z are for number density caculation;
#		 mass_xy calculates 2D mass density on XY plane;
#		 while number_x calculates 1D number density along X axis.
# firstFrame lastFrame step: specify the frame range and step
# xl xh xr: xl and xh define the two boundaries along X direction, xr is the number of slices between them;
#			so do the yl yh yr, and zl zh zr.
# atomSelection: user-defined selection, the text must be in {}
# length unit: angstrom; mass density unit g/cm^3; number density: /nm^3; total mass unit: amu


proc dc {file_1 file_2 type firstFrame lastFrame step \
         xl xh xr yl yh yr zl zh zr atomSelection} {

if {[string equal $file_1 $file_2]} {
	mol new $file_1 first $firstFrame last $lastFrame step $step waitfor all
} else {
	mol new $file_1
	mol addfile $file_2 first $firstFrame last $lastFrame step $step waitfor all
}

set dx [expr double($xh-$xl)/$xr]
set dy [expr double($yh-$yl)/$yr]
set dz [expr double($zh-$zl)/$zr]

# ------Mass density on YZ plane (unit:g/cm^3)------
if {[string equal -nocase $type mass_yz]} {
	for {set j 0} {$j<$xr} {incr j 1} {	
		set totalMass 0.0
		set filename [format "$type _%d_%d_%.3f_%.3f" $firstFrame $lastFrame [expr {$xl+$dx*$j}] [expr {$xl+$dx*($j+1)}]]
		set writefile [open $filename.dat w+]
		puts $writefile "------Mass density on YZ plane (unit:g/cm^3)------"
		puts $writefile [format "dx = %.3f, dy = %.3f, dz = %.3f" $dx $dy $dz]
		puts $writefile "xpos   ypos   zpos   massDensity"
		puts "------Mass density on YZ plane (unit:g/cm^3)------"
		puts [format "dx = %.3f, dy = %.3f, dz = %.3f" $dx $dy $dz]
		puts "xpos  ypos  zpos  massDensity"
		for {set k 0} {$k<$yr} {incr k 1} {
			for {set m 0} {$m<$zr} {incr m 1} {
				set bin [format "x>%.3f and x<%.3f and y>%.3f and y<%.3f and z>%.3f and z<%.3f" \
						[expr {$xl+$dx*$j}] [expr {$xl+$dx*($j+1)}] [expr {$yl+$dy*$k}] \
						[expr {$yl+$dy*($k+1)}] [expr {$zl+$dz*$m}] [expr {$zl+$dz*($m+1)}]]			
				set selection [atomselect top "$bin and $atomSelection"]
				set sum_frame 0
				set frameNub 0	
				for {set i $firstFrame} {$i<=$lastFrame} {incr i $step} {
					set frameNub [expr {$frameNub+1}]
					$selection frame $i
					$selection update
					foreach atomMass [$selection get mass] {
						set sum_frame [expr {$sum_frame+$atomMass}]
					}					
				}
				set xpos [expr {$xl+$dx*($j+0.5)}]
				set ypos [expr {$yl+$dy*($k+0.5)}]
				set zpos [expr {$zl+$dz*($m+0.5)}]
				set totalMass [expr {double($sum_frame)/$frameNub+$totalMass}]
				set massDensity [expr {double($sum_frame)*1.6605/($frameNub*$dx*$dy*$dz)}]
				puts $writefile [format "%.3f %.3f %.3f %.6f" $xpos $ypos $zpos $massDensity]
				puts [format "%.3f %.3f %.3f %.6f" $xpos $ypos $zpos $massDensity]	
			}
		}	
		puts $writefile [format "The total mass is %.6f amu" $totalMass]
		puts [format "The total mass is %.6f amu" $totalMass]
		close $writefile
	}
	
# ------Mass density on XZ plane (unit:g/cm^3)------	
} elseif {[string equal -nocase $type mass_xz]} {
	for {set j 0} {$j<$yr} {incr j 1} {	
		set totalMass 0.0
		set filename [format "$type _%d_%d_%.3f_%.3f" $firstFrame $lastFrame [expr {$yl+$dy*$j}] [expr {$yl+$dy*($j+1)}]]
		set writefile [open $filename.dat w+]
		puts $writefile "------Mass density on XZ plane (unit:g/cm^3)------"
		puts $writefile [format "dx = %.3f, dy = %.3f, dz = %.3f" $dx $dy $dz]
		puts $writefile "xpos   ypos   zpos   massDensity"
		puts "------Mass density on XZ plane (unit:g/cm^3)------"
		puts [format "dx = %.3f, dy = %.3f, dz = %.3f" $dx $dy $dz]
		puts "xpos  ypos  zpos  massDensity"
		for {set k 0} {$k<$xr} {incr k 1} {
			for {set m 0} {$m<$zr} {incr m 1} {
				set bin [format "x>%.3f and x<%.3f and y>%.3f and y<%.3f and z>%.3f and z<%.3f" \
						[expr {$xl+$dx*$k}] [expr {$xl+$dx*($k+1)}] [expr {$yl+$dy*$j}] \
						[expr {$yl+$dy*($j+1)}] [expr {$zl+$dz*$m}] [expr {$zl+$dz*($m+1)}]]			
				set selection [atomselect top "$bin and $atomSelection"]
				set sum_frame 0
				set frameNub 0	
				for {set i $firstFrame} {$i<=$lastFrame} {incr i $step} {
					set frameNub [expr {$frameNub+1}]
					$selection frame $i
					$selection update
					foreach atomMass [$selection get mass] {
						set sum_frame [expr {$sum_frame+$atomMass}]
					}					
				}
				set xpos [expr {$xl+$dx*($k+0.5)}]
				set ypos [expr {$yl+$dy*($j+0.5)}]
				set zpos [expr {$zl+$dz*($m+0.5)}]
				set totalMass [expr {double($sum_frame)/$frameNub+$totalMass}]
				set massDensity [expr {double($sum_frame)*1.6605/($frameNub*$dx*$dy*$dz)}]
				puts $writefile [format "%.3f %.3f %.3f %.6f" $xpos $ypos $zpos $massDensity]
				puts [format "%.3f %.3f %.3f %.6f" $xpos $ypos $zpos $massDensity]	
			}
		}	
		puts $writefile [format "The total mass is %.6f amu" $totalMass]
		puts [format "The total mass is %.6f amu" $totalMass]
		close $writefile
	}
	
# ------Mass density on XY plane (unit:g/cm^3)------	
} elseif {[string equal -nocase $type mass_xy]} {
	for {set j 0} {$j<$zr} {incr j 1} {	
		set totalMass 0.0
		set filename [format "$type _%d_%d_%.3f_%.3f" $firstFrame $lastFrame [expr {$zl+$dz*$j}] [expr {$zl+$dz*($j+1)}]]
		set writefile [open $filename.dat w+]
		puts $writefile "------Mass density on XY plane (unit:g/cm^3)------"
		puts $writefile [format "dx = %.3f, dy = %.3f, dz = %.3f" $dx $dy $dz]
		puts $writefile "xpos   ypos   zpos   massDensity"
		puts "------Mass density on XY plane (unit:g/cm^3)------"
		puts [format "dx = %.3f, dy = %.3f, dz = %.3f" $dx $dy $dz]
		puts "xpos  ypos  zpos  massDensity"
		for {set k 0} {$k<$xr} {incr k 1} {
			for {set m 0} {$m<$yr} {incr m 1} {
				set bin [format "x>%.3f and x<%.3f and y>%.3f and y<%.3f and z>%.3f and z<%.3f" \
						[expr {$xl+$dx*$k}] [expr {$xl+$dx*($k+1)}] [expr {$yl+$dy*$m}] \
						[expr {$yl+$dy*($m+1)}] [expr {$zl+$dz*$j}] [expr {$zl+$dz*($j+1)}]]			
				set selection [atomselect top "$bin and $atomSelection"]
				set sum_frame 0
				set frameNub 0	
				for {set i $firstFrame} {$i<=$lastFrame} {incr i $step} {
					set frameNub [expr {$frameNub+1}]
					$selection frame $i
					$selection update
					foreach atomMass [$selection get mass] {
						set sum_frame [expr {$sum_frame+$atomMass}]
					}					
				}
				set xpos [expr {$xl+$dx*($k+0.5)}]
				set ypos [expr {$yl+$dy*($m+0.5)}]
				set zpos [expr {$zl+$dz*($j+0.5)}]
				set totalMass [expr {double($sum_frame)/$frameNub+$totalMass}]
				set massDensity [expr {double($sum_frame)*1.6605/($frameNub*$dx*$dy*$dz)}]
				puts $writefile [format "%.3f %.3f %.3f %.6f" $xpos $ypos $zpos $massDensity]
				puts [format "%.3f %.3f %.3f %.6f" $xpos $ypos $zpos $massDensity]	
			}
		}	
		puts $writefile [format "The total mass is %.6f amu" $totalMass]
		puts [format "The total mass is %.6f amu" $totalMass]
		close $writefile
	}
	
# ------Mass density on X axis (unit:g/cm^3)------	
} elseif {[string equal -nocase $type mass_x]} {
	set dy [expr {$yh-$yl}]
	set dz [expr {$zh-$zl}]
	set totalMass 0.0
	set filename [format "$type _%d_%d" $firstFrame $lastFrame]
	set writefile [open $filename.dat w+]
	puts $writefile "------Mass density on X axis (unit:g/cm^3)------"
	puts $writefile [format "dx = %.3f, xl = %.3f, xh = %.3f, yl = %.3f, yh = %.3f, \
					 zl = %.3f, zh = %.3f" $dx $xl $xh $yl $yh $zl $zh]
	puts $writefile "xpos  massDensity"
	puts "------Mass density on X axis (unit:g/cm^3)------"
	puts [format "dx = %.3f, xl = %.3f, xh = %.3f, yl = %.3f, yh = %.3f, zl = %.3f, zh = %.3f"\
				 $dx $xl $xh $yl $yh $zl $zh]
	puts "xpos  massDensity"
	for {set j 0} {$j<$xr} {incr j 1} {
		set bin [format "x>%.3f and x<%.3f and y>%.3f and y<%.3f and z>%.3f and z<%.3f" \
				[expr {$xl+$dx*$j}] [expr {$xl+$dx*($j+1)}] $yl $yh $zl $zh]			
		set selection [atomselect top "$bin and $atomSelection"]
		set sum_frame 0
		set frameNub 0
		for {set i $firstFrame} {$i<=$lastFrame} {incr i $step} {
			set frameNub [expr {$frameNub+1}]
			$selection frame $i
			$selection update
			foreach atomMass [$selection get mass] {
				set sum_frame [expr {$sum_frame+$atomMass}]
			}					
		}
		set xpos [expr {$xl+$dx*($j+0.5)}]
		set totalMass [expr {double($sum_frame)/$frameNub+$totalMass}]
		set massDensity [expr {double($sum_frame)*1.6605/($frameNub*$dx*$dy*$dz)}]
		puts $writefile [format "%.3f %.6f" $xpos $massDensity]
		puts [format "%.3f %.6f" $xpos $massDensity]
	}
		puts $writefile [format "The total mass is %.6f amu" $totalMass]
		puts [format "The total mass is %.6f amu" $totalMass]
	close $writefile
	
# ------Mass density on Y axis (unit:g/cm^3)------	
} elseif {[string equal -nocase $type mass_y]} {
	set dx [expr {$xh-$xl}]
	set dz [expr {$zh-$zl}]
	set totalMass 0.0
	set filename [format "$type _%d_%d" $firstFrame $lastFrame]
	set writefile [open $filename.dat w+]
	puts $writefile "------Mass density on Y axis (unit:g/cm^3)------"
	puts $writefile [format "dy = %.3f, xl = %.3f, xh = %.3f, yl = %.3f, yh = %.3f, \
					 zl = %.3f, zh = %.3f" $dy $xl $xh $yl $yh $zl $zh]
	puts $writefile "ypos  massDensity"
	puts "------Mass density on Y axis (unit:g/cm^3)------"
	puts [format "dy = %.3f, xl = %.3f, xh = %.3f, yl = %.3f, yh = %.3f, zl = %.3f, zh = %.3f"\
				 $dy $xl $xh $yl $yh $zl $zh]
	puts "ypos  massDensity"
	for {set j 0} {$j<$yr} {incr j 1} {
		set bin [format "x>%.3f and x<%.3f and y>%.3f and y<%.3f and z>%.3f and z<%.3f" \
				$xl $xh [expr {$yl+$dy*$j}] [expr {$yl+$dy*($j+1)}]  $zl $zh]			
		set selection [atomselect top "$bin and $atomSelection"]
		set sum_frame 0
		set frameNub 0
		for {set i $firstFrame} {$i<=$lastFrame} {incr i $step} {
			set frameNub [expr {$frameNub+1}]
			$selection frame $i
			$selection update
			foreach atomMass [$selection get mass] {
				set sum_frame [expr {$sum_frame+$atomMass}]
			}					
		}
		set ypos [expr {$yl+$dy*($j+0.5)}]
		set totalMass [expr {double($sum_frame)/$frameNub+$totalMass}]
		set massDensity [expr {double($sum_frame)*1.6605/($frameNub*$dx*$dy*$dz)}]
		puts $writefile [format "%.3f %.6f" $ypos $massDensity]
		puts [format "%.3f %.6f" $ypos $massDensity]
	}
		puts $writefile [format "The total mass is %.6f amu" $totalMass]
		puts [format "The total mass is %.6f amu" $totalMass]
	close $writefile
	
# ------Mass density on Z axis (unit:g/cm^3)------	
} elseif {[string equal -nocase $type mass_z]} {
	set dy [expr {$yh-$yl}]
	set dx [expr {$xh-$xl}]
	set totalMass 0.0
	set filename [format "$type _%d_%d" $firstFrame $lastFrame]
	set writefile [open $filename.dat w+]
	puts $writefile "------Mass density on Z axis (unit:g/cm^3)------"
	puts $writefile [format "dz = %.3f, xl = %.3f, xh = %.3f, yl = %.3f, yh = %.3f, \
					 zl = %.3f, zh = %.3f" $dz $xl $xh $yl $yh $zl $zh]
	puts $writefile "zpos  massDensity"
	puts "------Mass density on Z axis (unit:g/cm^3)------"
	puts [format "dz = %.3f, xl = %.3f, xh = %.3f, yl = %.3f, yh = %.3f, zl = %.3f, zh = %.3f"\
				 $dz $xl $xh $yl $yh $zl $zh]
	puts "zpos  massDensity"
	for {set j 0} {$j<$zr} {incr j 1} {
		set bin [format "x>%.3f and x<%.3f and y>%.3f and y<%.3f and z>%.3f and z<%.3f" \
				$xl $xh $yl $yh [expr {$zl+$dz*$j}] [expr {$zl+$dz*($j+1)}]]			
		set selection [atomselect top "$bin and $atomSelection"]
		set sum_frame 0
		set frameNub 0
		for {set i $firstFrame} {$i<=$lastFrame} {incr i $step} {
			set frameNub [expr {$frameNub+1}]
			$selection frame $i
			$selection update
			foreach atomMass [$selection get mass] {
				set sum_frame [expr {$sum_frame+$atomMass}]
			}					
		}
		set zpos [expr {$zl+$dz*($j+0.5)}]
		set totalMass [expr {double($sum_frame)/$frameNub+$totalMass}]
		set massDensity [expr {double($sum_frame)*1.6605/($frameNub*$dx*$dy*$dz)}]
		puts $writefile [format "%.3f %.6f" $zpos $massDensity]
		puts [format "%.3f %.6f" $zpos $massDensity]
	}
		puts $writefile [format "The total mass is %.6f amu" $totalMass]
		puts [format "The total mass is %.6f amu" $totalMass]
	close $writefile
	
# ------Number density on YZ plane (unit:/nm^3)------	
} elseif {[string equal -nocase $type number_yz]} {
	for {set j 0} {$j<$xr} {incr j 1} {	
		set totalNumber 0
		set filename [format "$type _%d_%d_%.3f_%.3f" $firstFrame $lastFrame [expr {$xl+$dx*$j}] [expr {$xl+$dx*($j+1)}]]
		set writefile [open $filename.dat w+]
		puts $writefile "------Number density on YZ plane (unit:/nm^3)------"
		puts $writefile [format "dx = %.3f, dy = %.3f, dz = %.3f" $dx $dy $dz]
		puts $writefile "xpos   ypos   zpos   numberDensity"
		puts "------Number density on YZ plane (unit:/nm^3)------"
		puts [format "dx = %.3f, dy = %.3f, dz = %.3f" $dx $dy $dz]
		puts "xpos  ypos  zpos  numberDensity"
		for {set k 0} {$k<$yr} {incr k 1} {
			for {set m 0} {$m<$zr} {incr m 1} {
				set bin [format "x>%.3f and x<%.3f and y>%.3f and y<%.3f and z>%.3f and z<%.3f" \
						[expr {$xl+$dx*$j}] [expr {$xl+$dx*($j+1)}] [expr {$yl+$dy*$k}] \
						[expr {$yl+$dy*($k+1)}] [expr {$zl+$dz*$m}] [expr {$zl+$dz*($m+1)}]]			
				set selection [atomselect top "$bin and $atomSelection"]
				set sum_frame 0
				set frameNub 0	
				for {set i $firstFrame} {$i<=$lastFrame} {incr i $step} {
					set frameNub [expr {$frameNub+1}]
					$selection frame $i
					$selection update				
					set atomNumber [$selection num]
					set sum_frame [expr {$sum_frame+$atomNumber}]			
				}
				set xpos [expr {$xl+$dx*($j+0.5)}]
				set ypos [expr {$yl+$dy*($k+0.5)}]
				set zpos [expr {$zl+$dz*($m+0.5)}]
				set totalNumber [expr {$sum_frame/$frameNub+$totalNumber}]
				set numberDensity [expr {double($sum_frame)*1000/($frameNub*$dx*$dy*$dz)}]
				puts $writefile [format "%.3f %.3f %.3f %.6f" $xpos $ypos $zpos $numberDensity]
				puts [format "%.3f %.3f %.3f %.6f" $xpos $ypos $zpos $numberDensity]
			}
		}	
		puts $writefile [format "The total number is %d" $totalNumber]
		puts [format "The total number is %d" $totalNumber]
		close $writefile
	}
	
# ------Number density on XZ plane (unit:/nm^3)------
} elseif {[string equal -nocase $type number_xz]} {
	for {set j 0} {$j<$xr} {incr j 1} {	
		set totalNumber 0
		set filename [format "$type _%d_%d_%.3f_%.3f" $firstFrame $lastFrame [expr {$yl+$dy*$j}] [expr {$yl+$dy*($j+1)}]]
		set writefile [open $filename.dat w+]
		puts $writefile "------Number density on XZ plane (unit:/nm^3)------"
		puts $writefile [format "dx = %.3f, dy = %.3f, dz = %.3f" $dx $dy $dz]
		puts $writefile "xpos   ypos   zpos   numberDensity"
		puts "------Number density on XZ plane (unit:/nm^3)------"
		puts [format "dx = %.3f, dy = %.3f, dz = %.3f" $dx $dy $dz]
		puts "xpos  ypos  zpos  numberDensity"
		for {set k 0} {$k<$xr} {incr k 1} {
			for {set m 0} {$m<$zr} {incr m 1} {
				set bin [format "x>%.3f and x<%.3f and y>%.3f and y<%.3f and z>%.3f and z<%.3f" \
						[expr {$xl+$dx*$k}] [expr {$xl+$dx*($k+1)}] [expr {$yl+$dy*$j}] \
						[expr {$yl+$dy*($j+1)}] [expr {$zl+$dz*$m}] [expr {$zl+$dz*($m+1)}]]			
				set selection [atomselect top "$bin and $atomSelection"]
				set sum_frame 0
				set frameNub 0	
				for {set i $firstFrame} {$i<=$lastFrame} {incr i $step} {
					set frameNub [expr {$frameNub+1}]
					$selection frame $i
					$selection update				
					set atomNumber [$selection num]
					set sum_frame [expr {$sum_frame+$atomNumber}]			
				}
				set xpos [expr {$xl+$dx*($k+0.5)}]
				set ypos [expr {$yl+$dy*($j+0.5)}]
				set zpos [expr {$zl+$dz*($m+0.5)}]
				set totalNumber [expr {$sum_frame/$frameNub+$totalNumber}]
				set numberDensity [expr {double($sum_frame)*1000/($frameNub*$dx*$dy*$dz)}]
				puts $writefile [format "%.3f %.3f %.3f %.6f" $xpos $ypos $zpos $numberDensity]
				puts [format "%.3f %.3f %.3f %.6f" $xpos $ypos $zpos $numberDensity]
			}
		}	
		puts $writefile [format "The total number is %d" $totalNumber]
		puts [format "The total number is %d" $totalNumber]
		close $writefile
	}
	
# ------Number density on XY plane (unit:/nm^3)------
} elseif {[string equal -nocase $type number_xy]} {
	for {set j 0} {$j<$zr} {incr j 1} {	
		set totalNumber 0
		set filename [format "$type _%d_%d_%.3f_%.3f" $firstFrame $lastFrame [expr {$zl+$dz*$j}] [expr {$zl+$dz*($j+1)}]]
		set writefile [open $filename.dat w+]
		puts $writefile "------Number density on XY plane (unit:/nm^3)------"
		puts $writefile [format "dx = %.3f, dy = %.3f, dz = %.3f" $dx $dy $dz]
		puts $writefile "xpos   ypos   zpos   numberDensity"
		puts "------Number density on XY plane (unit:/nm^3)------"
		puts [format "dx = %.3f, dy = %.3f, dz = %.3f" $dx $dy $dz]
		puts "xpos  ypos  zpos  numberDensity"
		for {set k 0} {$k<$yr} {incr k 1} {
			for {set m 0} {$m<$xr} {incr m 1} {
				set bin [format "x>%.3f and x<%.3f and y>%.3f and y<%.3f and z>%.3f and z<%.3f" \
						[expr {$xl+$dx*$m}] [expr {$xl+$dx*($m+1)}] [expr {$yl+$dy*$k}] \
						[expr {$yl+$dy*($k+1)}] [expr {$zl+$dz*$j}] [expr {$zl+$dz*($j+1)}]]			
				set selection [atomselect top "$bin and $atomSelection"]
				set sum_frame 0
				set frameNub 0	
				for {set i $firstFrame} {$i<=$lastFrame} {incr i $step} {
					set frameNub [expr {$frameNub+1}]
					$selection frame $i
					$selection update				
					set atomNumber [$selection num]
					set sum_frame [expr {$sum_frame+$atomNumber}]			
				}
				set xpos [expr {$xl+$dx*($m+0.5)}]
				set ypos [expr {$yl+$dy*($k+0.5)}]
				set zpos [expr {$zl+$dz*($j+0.5)}]
				set totalNumber [expr {$sum_frame/$frameNub+$totalNumber}]
				set numberDensity [expr {double($sum_frame)*1000/($frameNub*$dx*$dy*$dz)}]
				puts $writefile [format "%.3f %.3f %.3f %.6f" $xpos $ypos $zpos $numberDensity]
				puts [format "%.3f %.3f %.3f %.6f" $xpos $ypos $zpos $numberDensity]
			}
		}	
		puts $writefile [format "The total number is %d" $totalNumber]
		puts [format "The total number is %d" $totalNumber]
		close $writefile
	}
	
# ------Number density on X axis (unit:/nm^3)------
} elseif {[string equal -nocase $type number_x]} {
	set dy [expr {$yh-$yl}]
	set dz [expr {$zh-$zl}]
	set totalNumber 0
	set filename [format "$type _%d_%d" $firstFrame $lastFrame]
	set writefile [open $filename.dat w+]
	puts $writefile "------Number density on X axis (unit:/nm^3)------"
	puts $writefile [format "dx = %.3f, xl = %.3f, xh = %.3f, yl = %.3f, yh = %.3f, \
					 zl = %.3f, zh = %.3f" $dx $xl $xh $yl $yh $zl $zh]
	puts $writefile "xpos  numberDensity"
	puts "------Number density on X axis (unit:/nm^3)------"
	puts [format "dx = %.3f, xl = %.3f, xh = %.3f, yl = %.3f, yh = %.3f, zl = %.3f, zh = %.3f"\
				 $dx $xl $xh $yl $yh $zl $zh]
	puts "xpos  numberDensity"
	for {set j 0} {$j<$xr} {incr j 1} {
		set bin [format "x>%.3f and x<%.3f and y>%.3f and y<%.3f and z>%.3f and z<%.3f" \
				[expr {$xl+$dx*$j}] [expr {$xl+$dx*($j+1)}] $yl $yh $zl $zh]			
		set selection [atomselect top "$bin and $atomSelection"]
		set sum_frame 0
		set frameNub 0
		for {set i $firstFrame} {$i<=$lastFrame} {incr i $step} {
			set frameNub [expr {$frameNub+1}]
			$selection frame $i
			$selection update				
			set atomNumber [$selection num]
			set sum_frame [expr {$sum_frame+$atomNumber}]			
		}
		set xpos [expr {$xl+$dx*($j+0.5)}]
		set totalNumber [expr {$sum_frame/$frameNub+$totalNumber}]
		set numberDensity [expr {double($sum_frame)*1000/($frameNub*$dx*$dy*$dz)}]
		puts $writefile [format "%.3f %.6f" $xpos $numberDensity]
		puts [format "%.3f %.6f" $xpos $numberDensity]
	}
		puts $writefile [format "The total number is %d" $totalNumber]
		puts [format "The total number is %d" $totalNumber]
	close $writefile
	
# ------Number density on Y axis (unit:/nm^3)------
} elseif {[string equal -nocase $type number_y]} {
	set dx [expr {$xh-$xl}]
	set dz [expr {$zh-$zl}]
	set totalNumber 0
	set filename [format "$type _%d_%d" $firstFrame $lastFrame]
	set writefile [open $filename.dat w+]
	puts $writefile "------Number density on Y axis (unit:/nm^3)------"
	puts $writefile [format "dy = %.3f, xl = %.3f, xh = %.3f, yl = %.3f, yh = %.3f, \
					 zl = %.3f, zh = %.3f" $dy $xl $xh $yl $yh $zl $zh]
	puts $writefile "ypos  numberDensity"
	puts "------Number density on Y axis (unit:/nm^3)------"
	puts [format "dy = %.3f, xl = %.3f, xh = %.3f, yl = %.3f, yh = %.3f, zl = %.3f, zh = %.3f"\
				 $dy $xl $xh $yl $yh $zl $zh]
	puts "ypos  numberDensity"
	for {set j 0} {$j<$yr} {incr j 1} {
		set bin [format "x>%.3f and x<%.3f and y>%.3f and y<%.3f and z>%.3f and z<%.3f" \
				$xl $xh [expr {$yl+$dy*$j}] [expr {$yl+$dy*($j+1)}]  $zl $zh]			
		set selection [atomselect top "$bin and $atomSelection"]
		set sum_frame 0
		set frameNub 0
		for {set i $firstFrame} {$i<=$lastFrame} {incr i $step} {
			set frameNub [expr {$frameNub+1}]
			$selection frame $i
			$selection update
			set atomNumber [$selection num]
			set sum_frame [expr {$sum_frame+$atomNumber}]			
		}
		set ypos [expr {$yl+$dy*($j+0.5)}]
		set totalNumber [expr {$sum_frame/$frameNub+$totalNumber}]
		set numberDensity [expr {double($sum_frame)*1000/($frameNub*$dx*$dy*$dz)}]
		puts $writefile [format "%.3f %.6f" $ypos $numberDensity]
		puts [format "%.3f %.6f" $ypos $numberDensity]
	}
		puts $writefile [format "The total number is %d" $totalNumber]
		puts [format "The total number is %d" $totalNumber]
	close $writefile
	
# ------Number density on Z axis (unit:/nm^3)------
}  elseif {[string equal -nocase $type number_z]} {
	set dy [expr {$yh-$yl}]
	set dx [expr {$xh-$xl}]
	set totalNumber 0
	set filename [format "$type _%d_%d" $firstFrame $lastFrame]
	set writefile [open $filename.dat w+]
	puts $writefile "------Number density on Z axis (unit:/nm^3)------"
	puts $writefile [format "dz = %.3f, xl = %.3f, xh = %.3f, yl = %.3f, yh = %.3f, \
					 zl = %.3f, zh = %.3f" $dz $xl $xh $yl $yh $zl $zh]
	puts $writefile "zpos  numberDensity"
	puts "------Number density on Z axis (unit:/nm^3)------"
	puts [format "dz = %.3f, xl = %.3f, xh = %.3f, yl = %.3f, yh = %.3f, zl = %.3f, zh = %.3f"\
				 $dz $xl $xh $yl $yh $zl $zh]
	puts "zpos  numberDensity"
	for {set j 0} {$j<$zr} {incr j 1} {
		set bin [format "x>%.3f and x<%.3f and y>%.3f and y<%.3f and z>%.3f and z<%.3f" \
				$xl $xh $yl $yh [expr {$zl+$dz*$j}] [expr {$zl+$dz*($j+1)}] ]			
		set selection [atomselect top "$bin and $atomSelection"]
		set sum_frame 0
		set frameNub 0
		for {set i $firstFrame} {$i<=$lastFrame} {incr i $step} {
			set frameNub [expr {$frameNub+1}]
			$selection frame $i
			$selection update				
			set atomNumber [$selection num]
			set sum_frame [expr {$sum_frame+$atomNumber}]			
		}
		set zpos [expr {$zl+$dz*($j+0.5)}]
		set totalNumber [expr {$sum_frame/$frameNub+$totalNumber}]
		set numberDensity [expr {double($sum_frame)*1000/($frameNub*$dx*$dy*$dz)}]
		puts $writefile [format "%.3f %.6f" $zpos $numberDensity]
		puts [format "%.3f %.6f" $zpos $numberDensity]
	}
		puts $writefile [format "The total number is %d" $totalNumber]
		puts [format "The total number is %d" $totalNumber]
	close $writefile
}

}  
# end of proc