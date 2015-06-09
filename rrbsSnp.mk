##
## Auto Generated makefile by CodeLite IDE
## any manual changes will be erased      
##
## Release
ProjectName            :=rrbsSnp
IntermediateDirectory  :=./
LinkerName             :=g++
ObjectSuffix           :=.o
DependSuffix           :=.o.d
PreprocessSuffix       :=.i
IncludeSwitch          :=-I
LibrarySwitch          :=-l
OutputSwitch           :=-o 
LibraryPathSwitch      :=-L
PreprocessorSwitch     :=-D
SourceSwitch           :=-c 
OutputFile             :=./$(ProjectName)
Preprocessors          :=
ObjectSwitch           :=-o 
ArchiveOutputSwitch    := 
PreprocessOnlySwitch   :=-E
ObjectsFileList        :="rrbsSnp.txt"
PCHCompileFlags        :=
MakeDirCommand         :=mkdir
RcCmpOptions           := 
LinkOptions            :=  
IncludePath            :=$(IncludeSwitch). $(IncludeSwitch). 
IncludePCH             := 
RcIncludePath          := 
Libs                   := 
ArLibs                 :=  
LibPath                :=$(LibraryPathSwitch). -L./samtools-0.1.19/ -lbam -lz -lpthread

##
## Common variables
## AR, CXX, CC, AS, CXXFLAGS and CFLAGS can be overriden using an environment variables
##
AR       := ar rcu
CXX      := g++
CC       := g++ 
CXXFLAGS := -O2 -static -m64 -I./samtools-0.1.19/ -L./samtools-0.1.19/ -lbam -lz -lpthread
CFLAGS   := -O2 -static -m64 -I./samtools-0.1.19/ -L./samtools-0.1.19/ -lbam -lz -lpthread
ASFLAGS  := 
AS       := as
RM		 := rm

##
## User defined environment variables
##
Objects0=$(IntermediateDirectory)/main.c$(ObjectSuffix) $(IntermediateDirectory)/chrome_funcs.c$(ObjectSuffix) $(IntermediateDirectory)/hash_funcs.c$(ObjectSuffix) $(IntermediateDirectory)/sam_funcs.c$(ObjectSuffix) 

Objects=$(Objects0) 

##
## Main Build Targets 
##
.PHONY: all clean PreBuild PrePreBuild PostBuild
all: $(OutputFile)

$(OutputFile): $(Objects) 
	@echo $(Objects0)  > $(ObjectsFileList)
	$(LinkerName) $(OutputSwitch)$(OutputFile) @$(ObjectsFileList) $(LibPath) $(Libs) $(LinkOptions)

PreBuild:


##
## Objects
##
$(IntermediateDirectory)/main.c$(ObjectSuffix): main.c $(IntermediateDirectory)/main.c$(DependSuffix)
	$(CC) $(SourceSwitch) "./main.c" $(CFLAGS) $(ObjectSwitch)$(IntermediateDirectory)/main.c$(ObjectSuffix) $(IncludePath)
$(IntermediateDirectory)/main.c$(DependSuffix): main.c
	@$(CC) $(CFLAGS) $(IncludePath) -MG -MP -MT$(IntermediateDirectory)/main.c$(ObjectSuffix) -MF$(IntermediateDirectory)/main.c$(DependSuffix) -MM "main.c"

$(IntermediateDirectory)/main.c$(PreprocessSuffix): main.c
	@$(CC) $(CFLAGS) $(IncludePath) $(PreprocessOnlySwitch) $(OutputSwitch) $(IntermediateDirectory)/main.c$(PreprocessSuffix) "main.c"

$(IntermediateDirectory)/chrome_funcs.c$(ObjectSuffix): chrome_funcs.c $(IntermediateDirectory)/chrome_funcs.c$(DependSuffix)
	$(CC) $(SourceSwitch) "./chrome_funcs.c" $(CFLAGS) $(ObjectSwitch)$(IntermediateDirectory)/chrome_funcs.c$(ObjectSuffix) $(IncludePath)
$(IntermediateDirectory)/chrome_funcs.c$(DependSuffix): chrome_funcs.c
	@$(CC) $(CFLAGS) $(IncludePath) -MG -MP -MT$(IntermediateDirectory)/chrome_funcs.c$(ObjectSuffix) -MF$(IntermediateDirectory)/chrome_funcs.c$(DependSuffix) -MM "chrome_funcs.c"

$(IntermediateDirectory)/chrome_funcs.c$(PreprocessSuffix): chrome_funcs.c
	@$(CC) $(CFLAGS) $(IncludePath) $(PreprocessOnlySwitch) $(OutputSwitch) $(IntermediateDirectory)/chrome_funcs.c$(PreprocessSuffix) "chrome_funcs.c"

$(IntermediateDirectory)/hash_funcs.c$(ObjectSuffix): hash_funcs.c $(IntermediateDirectory)/hash_funcs.c$(DependSuffix)
	$(CC) $(SourceSwitch) "./hash_funcs.c" $(CFLAGS) $(ObjectSwitch)$(IntermediateDirectory)/hash_funcs.c$(ObjectSuffix) $(IncludePath)
$(IntermediateDirectory)/hash_funcs.c$(DependSuffix): hash_funcs.c
	@$(CC) $(CFLAGS) $(IncludePath) -MG -MP -MT$(IntermediateDirectory)/hash_funcs.c$(ObjectSuffix) -MF$(IntermediateDirectory)/hash_funcs.c$(DependSuffix) -MM "hash_funcs.c"

$(IntermediateDirectory)/hash_funcs.c$(PreprocessSuffix): hash_funcs.c
	@$(CC) $(CFLAGS) $(IncludePath) $(PreprocessOnlySwitch) $(OutputSwitch) $(IntermediateDirectory)/hash_funcs.c$(PreprocessSuffix) "hash_funcs.c"

$(IntermediateDirectory)/sam_funcs.c$(ObjectSuffix): sam_funcs.c $(IntermediateDirectory)/sam_funcs.c$(DependSuffix)
	$(CC) $(SourceSwitch) "./sam_funcs.c" $(CFLAGS) $(ObjectSwitch)$(IntermediateDirectory)/sam_funcs.c$(ObjectSuffix) $(IncludePath)
$(IntermediateDirectory)/sam_funcs.c$(DependSuffix): sam_funcs.c
	@$(CC) $(CFLAGS) $(IncludePath) -MG -MP -MT$(IntermediateDirectory)/sam_funcs.c$(ObjectSuffix) -MF$(IntermediateDirectory)/sam_funcs.c$(DependSuffix) -MM "sam_funcs.c"

$(IntermediateDirectory)/sam_funcs.c$(PreprocessSuffix): sam_funcs.c
	@$(CC) $(CFLAGS) $(IncludePath) $(PreprocessOnlySwitch) $(OutputSwitch) $(IntermediateDirectory)/sam_funcs.c$(PreprocessSuffix) "sam_funcs.c"


-include $(IntermediateDirectory)/*$(DependSuffix)
##
## Clean
##
clean:
	$(RM) ./*$(ObjectSuffix)
	$(RM) ./*$(DependSuffix)
	$(RM) $(OutputFile)




