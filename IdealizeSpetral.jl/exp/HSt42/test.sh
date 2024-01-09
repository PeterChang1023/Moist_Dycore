#!bin/sh


start_day=0
final_day=1000
space_day=100

for i in `seq $start_day $space_day $final_day`
do
L=25
echo $i"day"
echo -n $L > Latent_heat.txt

if [ $i -eq 0 ]; then
	rm -rf HSt42_${L}
	mkdir HSt42_${L}
	echo -n $space_day > HSt42_${L}/day_interval.txt
	echo -n "None" > HSt42_${L}/firstday_file.txt
else
	echo -n "warmstart.dat" > HSt42_${L}/firstday_file.txt
fi


julia Run_HS.jl

if [ -f "HSt42_${L}/warmstart.dat" ] && [ $i -lt $final_day ]; then
	cp "HSt42_"${L}"/warmstart.dat" "warmstart.dat"
	mv "HSt42_"${L}"/all_L"${L}".dat" "HSt42_"${L}"/RH80_PR"$L"_"$final_day"day_startfrom_"$i"day_final.dat"
	echo 'warmstart file exists.'
elif [ -f "HSt42_${L}/warmstart.dat" ] && [ $i -eq $final_day ]; then
	mv "HSt42_"${L}"/all_L"${L}".dat" "HSt42_"${L}"/RH80_PR"$L"_"$final_day"day_startfrom_"$i"day_final.dat"
	mv "HSt42_"${L}"/warmstart.dat" "HSt42_"${L}"/HSt42_"${L}"RH80_PR"$L"_"$final_day"day_startfrom_"$i"day_final.dat"
	echo "All files have completed!!!"
fi
done
