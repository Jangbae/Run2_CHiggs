#!/bin/bash
# 
echo "HADDING NOMINAL"
 
python -u haddoutput.py nominal >& haddnominal.log &

# echo "HADDING JECUP"
# 
# python -u haddoutput.py JECup >& haddJECup.log &
# 
# echo "HADDING JECDOWN"
# python -u haddoutput.py JECdown >& haddJECdown.log &
# 
# echo "HADDING JERUP"
# 
# python -u haddoutput.py JERup >& haddJERup.log &
# 
# echo "HADDING JERDOWN"
# 
# python -u haddoutput.py JERdown >& haddJERdown.log &

echo "DONE"
