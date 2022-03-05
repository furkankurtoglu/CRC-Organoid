import smtplib, ssl

smtp_server = "smtp.gmail.com"
port = 587  # For starttls
sender_email = "pcnotifier20@gmail.com"
password = "pc2020iu"

# Create a secure SSL context
context = ssl.create_default_context()

# Try to log in to server and send email
try:
    server = smtplib.SMTP(smtp_server,port)
    server.ehlo() # Can be omitted
    server.starttls(context=context) # Secure the connection
    server.ehlo() # Can be omitted
    server.login(sender_email, password)
    sender_email = "pcnotifier20@gmail.com"
    receiver_email = "f.kurtoglu9@gmail.com"
    message = "Subject: Run Completed! \nDear PhysiCell user,\n\nYour simulation has been completed succesfully.\nHave a nice day! \n\n---PhysiCell Notifier---"
    server.sendmail(sender_email, receiver_email, message)
    
except Exception as e:
    # Print any error messages to stdout
    print(e)
finally:
    server.quit()