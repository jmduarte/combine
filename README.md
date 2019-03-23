HiggsAnalysis-CombinedLimit
===========================

### Official documentation

[Manual to run combine](https://github.com/cms-analysis/HiggsAnalysis-CombinedLimit/wiki)

### Standalone compilation in `lxplus`
```
git clone https://github.com/cms-analysis/HiggsAnalysis-CombinedLimit.git combine
cd combine
source env_standalone.sh 
make -j 8; make # second make fixes compilation error of first
```
### Code formatting

To format the whole C++ code:
```
find . -regex '.*\.\(cpp\|hpp\|cc\|cxx\|hh\|hxx\|h\)' -exec clang-format -style=file -i {} \;
```
