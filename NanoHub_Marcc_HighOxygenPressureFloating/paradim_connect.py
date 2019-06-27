import paramiko
import logging
import socket
import time, threading
import re

class ShellHandler:

    def __init__(self, val_code, val_user, val_pwd, val_server, val_port):        
        def inter_handler(title, instructions, prompt_list):
            resp = [] 
            for pr in prompt_list:
                if str(pr[0]).strip() == "Username:":
                    resp.append(val_user)
                elif str(pr[0]).strip() == "Password:":
                    resp.append(val_pwd)
                elif str(pr[0]).strip() == "Verification code:":
                    resp.append(val_code)
            return tuple(resp)

        self.sock = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
        self.sock.connect((val_server, val_port))
        ts = paramiko.Transport(self.sock)
        ts.connect()
        ts.auth_interactive_dumb(val_user, inter_handler)
        self.ssh = ts.open_session()
        self.ssh.invoke_shell()
        self.stdin = self.ssh.makefile('wb')
        self.stdout = self.ssh.makefile('r')  
        self.stderr = self.ssh.makefile_stderr('r')  
        self.initialized = False       

    def __del__(self):
        try:
            self.ssh.close()
            #self.sock.close()
        except:
            pass

    def is_active(self):
        return not self.sock._closed

    def execute(self, cmd, buffer=None):
        cmd = cmd.strip('\n')
        #print (cmd)
        self.stdin.write(cmd + '\n')
        finish = 'end of stdOUT buffer. finished with exit status'
        echo_cmd = 'echo {} $?'.format(finish)
        self.stdin.write(echo_cmd + '\n')
        shin = self.stdin
        self.stdin.flush()

        shout = []
        sherr = []
        exit_status = 0
        for line in self.stdout:
            if str(line).startswith(cmd) or str(line).startswith(echo_cmd):
                # up for now filled with shell junk from stdin
                shout = []
            elif str(line).startswith(finish):
                # our finish command ends with the exit status
                exit_status = int(str(line).rsplit(maxsplit=1)[1])
                if exit_status:
                    # stderr is combined with stdout.
                    # thus, swap sherr with shout in a case of failure.
                    sherr = shout
                    shout = []
                break
            else:
                # get rid of 'coloring and formatting' special characters
                shout.append(re.compile(r'(\x9B|\x1B\[)[0-?]*[ -/]*[@-~]').sub('', line).
                             replace('\b', '').replace('\r', ''))

        # first and last lines of shout/sherr contain a prompt
        if shout and echo_cmd in shout[-1]:
            shout.pop()
        if shout and cmd in shout[0]:
            shout.pop(0)
        if sherr and echo_cmd in sherr[-1]:
            sherr.pop()
        if sherr and cmd in sherr[0]:
            sherr.pop(0)
            
        if buffer is not None:
            buffer['commands'].value = buffer['commands'].value + cmd + '\n'
            buffer['stdin'].value = buffer['stdin'].value + cmd + '\n'
            buffer['stdout'].value = buffer['stdout'].value + ''.join(shout)
            buffer['stderr'].value = buffer['stderr'].value + ''.join(sherr)
        return shin, shout, sherr, cmd


class setInterval :
    def __init__(self,interval,action) :
        self.interval=interval
        self.action=action
        self.stopEvent=threading.Event()
        thread=threading.Thread(target=self.__setInterval)
        thread.start()

    def __setInterval(self) :
        nextTime=time.time()+self.interval
        while not self.stopEvent.wait(nextTime-time.time()) :
            nextTime+=self.interval
            self.action()

    def cancel(self) :
        self.stopEvent.set()
