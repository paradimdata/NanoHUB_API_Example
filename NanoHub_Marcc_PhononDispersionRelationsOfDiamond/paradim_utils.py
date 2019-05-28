def GetPbs( template, input_file, command='pw.x' ):
    change = {  'reservation':'paradim', 
                'job_name':input_file,
                'time':'00:30:00',
                'partition':'shared', #debug
                'nodes':'1',
                'tasks':'12',
                'mem':'1000MB',
                'modules':'ml intel intelmpi quantumespresso/6.4.1',
                'command':'mpiexec -np 12 '+ command +' -npool 4 < ' + input_file + '.in >' + input_file + '.out '}
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

def isfloat(value):
  try:
    float(value)
    return True
  except ValueError:
    return False