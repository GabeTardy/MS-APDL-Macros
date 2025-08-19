@rem enough. just going to delete all the lock files.
del /q "*.lock"

@rem run ansys
ansys252 -g -dir %~dp0 -d win32c -smp