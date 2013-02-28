#!/bin/bash

# Modify these if you want.
num_threads=11

# These things should stay the same.
num_cities=17
program_name=wsp_openmp
num_iterations=3

# The head node has hostnames like tg-login1.blacklight.psc.teragrid.org,
# the batch job has hostnames like bl1.psc.teragrid.org. Searching for
# blacklight in the hostname is a poor mans filter for students running on the
# head node.
hostname | grep -v blacklight >/dev/null
if [ $? -ne 0 ]
then
echo "Do not run on blacklight head node!"
exit -1
fi

# Make sure the code is compiled and runs before testing.
make run
if [ $? -ne 0 ]
then
exit $?
fi

# Make sure they already downloaded a scoreboard token.
if [ ! -f ./scoreboard_token ]
then
  cat <<EOM
You must first download a scoreboard token to check your code.
Get one from http://dolores.sp.cs.cmu.edu/15418_spr13/index.php/scoreboard/token
Enter it here or save it to a file named "scoreboard_token".
EOM
  echo -n ": "
  read token
  if [ ${#token} -le 40 ]
  then
    echo "Invalid token: did you include your username?"
    exit 1
  else
    echo $token > scoreboard_token
  fi
fi

for seed in 15418 1000000 3735928559 111111111 11
do
  echo "Checking seed $seed..."
  # This should be big enough, right?
  best_time=9999999.0
  solution=-1
  infile=input/dist$num_cities.$seed
  python mkinput.py $num_cities --seed=$seed > $infile

  # If mkinput.py fails, we should fail.
  if [ $? -ne 0 ]
  then
    exit $?
  fi

  # Find the min of a number of iterations to rule out noise.
  for iter in `seq $num_iterations`
  do
    result=$(OMP_NUM_THREADS=$num_threads ./wsp -i $infile)
    # If ./wsp fails, we should fail.
    if [ $? -ne 0 ]
    then
      exit $?
    fi

    # This was really the point where I should have switched to python,
    # but gosh-darned-it I'm stubborn. The awk crap is necessary because
    # floating point is hard in bash.
    best_time=$(echo "$result" | grep took                       \
                                 | tr -d s                         \
                                 | awk "{                          \
                                           if(\$3 < $best_time) {  \
                                              print \$3            \
                                           } else {                \
                                               print $best_time    \
                                           }                       \
                                        }" 2>/dev/null)

    # We could check the solution every time, but we're lazy and will only
    # report the last one.
    solution=$(echo "$result" | grep distance | sed 's/.* of distance //')

    echo "After $iter iterations, best time for $seed is $best_time seconds"
  done
  echo "Sending best time to server..."
  curl -F token=$(cat ./scoreboard_token) \
       -F program_name=$program_name      \
       -F score=$best_time                \
       -F machine=$(hostname)             \
       -F instance=$seed                  \
       -F cores=$num_threads              \
       -F solution=$solution   \
       http://dolores.sp.cs.cmu.edu/15418_spr13/index.php/scoreboard/submit
done
