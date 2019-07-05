'''def GetPbs( template, input_file, command='pw.x' ):
    change = {  'reservation':'paradim', 
                'job_name':input_file,
                'time':'00:30:00',
                'partition':'shared', #debug
                'nodes':'1',
                'tasks':'24',
                'mem':'1000MB',
                'modules':'ml intel intelmpi quantumespresso/6.4.1',
                'command':'mpiexec -np 12 '+ command +' -npool 4 < ' + input_file + '.in >' + input_file + '.out '}
    return template.safe_substitute(change)
''' 

def GetPbs( template, input_file, command ):
    change = {  'reservation':'paradim', 
                'job_name':input_file,
                'time':'00:30:00',
                'partition':'shared', #debug
                'nodes':'1',
                'tasks':'12',
                'mem':'1000MB',
                'modules':'ml intel intelmpi quantumespresso/6.4.1',
                'command':command}
    return template.safe_substitute(change)

def getStatus( session, code, buffer = None ):
    command = "qstat -a " + str(code)
    try :
        stdin, stdout, stderr, command = session.execute(command, buffer);
        line = stdout[len(stdout)-1].strip('\n')
        columns = line.split()
        if (len(columns) < 11):
            status = 'X'
        else :
            status = str(columns[9])
    except :
        status = 'X'
        line = "ERROR"
    return status, line, command

def fileExist( session, file ):
    try :
        stdin, stdout, stderr, command = session.execute("test -e " + file + " && echo file exists || echo file not found");
        stdout = "".join(stdout).strip('\n')
        if (stdout != "file exists"):
            return False
        return True
    except :
        return False

def getUPF( session, file ):    
    if (fileExist(session, file) is False):
        session.execute("wget http://www.quantum-espresso.org/upf_files/" + v["upf"], PARADIM_UI['s1']);

    
def isfloat(value):
  try:
    float(value)
    return True
  except ValueError:
    return False