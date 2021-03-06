#include "gtest/gtest.h"
#include <unistd.h>
#include <fstream>
#include <string>
#include <vector>
#include <IRA>

#include "Schedule.h"

using namespace Schedule;

class TestProcedure : public ::testing::Test
{
    public:
        TestProcedure();
        virtual ~TestProcedure();
        IRA::CString filename, basepath;
        std::string fullname;
        char* _basepath;
        CProcedureList *cpl;
};

TestProcedure::TestProcedure() : 
    filename("procedure.cfg")
{
    char* _basepath = get_current_dir_name();
    basepath = IRA::CString(_basepath);
    fullname = std::string((const char*)(basepath + "/" + filename));
    cpl = new CProcedureList(basepath + "/", filename);
    free(_basepath);
}

TestProcedure::~TestProcedure()
{
    delete cpl;
}

TEST_F(TestProcedure, can_open_procedure_file)
{
    std::ifstream input_file;
    input_file.open(fullname.c_str());
    ASSERT_TRUE(input_file.good());
    ASSERT_TRUE(input_file.is_open());
    input_file.close();
}

TEST_F(TestProcedure, CProcedureList_constructor)
{
    ASSERT_EQ(cpl->getTotalLines(), (DWORD)0);
}

TEST_F(TestProcedure, CProcedureList_readAll)
{
    ASSERT_TRUE(cpl->readAll(false))
        << "error: " << (const char*)cpl->getLastError();
}

TEST_F(TestProcedure, CProcedureList_readAll_with_check)
{
    ASSERT_TRUE(cpl->readAll(true))
        << "error: " << (const char*)cpl->getLastError();
}

TEST_F(TestProcedure, CProcedureList_getProcedure_by_name)
{
    cpl->readAll(false);
    std::vector<IRA::CString> commands;
    WORD args;
    bool result = cpl->getProcedure("PROCEDURE_WAIT2", commands, args);
    ASSERT_TRUE(result) 
        << "error: " << (const char*)cpl->getLastError();
    ASSERT_EQ(commands[0], IRA::CString("wait=2"));
}

TEST_F(TestProcedure, CProcedureList_getProcedure_by_name_with_params)
{
    cpl->readAll(false);
    std::vector<IRA::CString> commands;
    WORD args;
    bool result = cpl->getProcedure("PROCEDURE_WAIT_PARAM", commands, args);
    ASSERT_TRUE(result) 
        << "error: " << (const char*)cpl->getLastError();
    ASSERT_EQ(args, (WORD)1) 
        << "error: " << (const char*)cpl->getLastError();
    ASSERT_EQ(commands[0], IRA::CString("wait=$1"));
}

TEST_F(TestProcedure,
CProcedureList_getProcedure_by_name_with_params_to_ACSStringSeq)
{
    cpl->readAll(false);
    ACS::stringSeq  commands;
    WORD args;
    bool result = cpl->getProcedure("PROCEDURE_WAIT_PARAM", commands, args);
    ASSERT_TRUE(result) 
        << "error: " << (const char*)cpl->getLastError();
    ASSERT_EQ(args, (WORD)1) 
        << "error: " << (const char*)cpl->getLastError();
    IRA::CString line((const char*)commands[0]);
    ASSERT_EQ(line, IRA::CString("wait=$1"));
}

TEST_F(TestProcedure, replace_procedure_parameters)
{
    const char PARAM = '$';
    cpl->readAll(false);
    ACS::stringSeq  commands;
    WORD args;
    cpl->getProcedure("PROCEDURE_WAIT_PARAM", commands, args);
    IRA::CString line((const char*)commands[0]);
    /* This acts as in Parser Lib, this test should be moved there for
     * consistency */
    IRA::CString params[1] = { IRA::CString("10") };
    IRA::CString param_string;
    param_string.Format("%c%d", PARAM, 1);
    line.ReplaceAll((const char*)param_string,(const char*)params[0]); 
    ASSERT_EQ(line, IRA::CString("wait=10"));
}

