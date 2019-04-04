from galaxy.jobs import JobDestination
import os

def setSGEpriority(app, user, user_email):
    admin_users = app.config.get( "admin_users", "" ).split( "," )

    try:
        roles = [role.name for role in user.all_roles()]
    except:
        roles = []

    params = {}
    
    """
    User priorities from TOP to BOTTOM.
    """
    # ADMINS
    if user_email in admin_users:
        params["nativeSpecification"] = "-p 1000"
    
    # TEACHERS
    elif "teacher" in roles:
        params["nativeSpecification"] = "-p 500"
    
    # PUPILS
    elif "pupil" in roles:
        params["nativeSpecification"] = "-p 10"
    
    # REGISTERED USERS
    elif user_email:
        params["nativeSpecification"] = "-p -100"
    
    # UNREGISTERED USERS
    else:
        params["nativeSpecification"] = "-p -600"


    return JobDestination(runner="drmaa", params=params)
