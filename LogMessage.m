function LogMessage(app, message)
    app.LogoutputTextArea.Value = [app.LogoutputTextArea.Value; {message}];
    drawnow;
end