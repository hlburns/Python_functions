
# coding: utf-8

# In[1]:

from IPython.html import widgets # Loads the Widget framework.
from IPython.core.magics.namespace import NamespaceMagics # Used to query namespace.

# For this example, hide these names, just to avoid polluting the namespace further
get_ipython().user_ns_hidden['widgets'] = widgets
get_ipython().user_ns_hidden['NamespaceMagics'] = NamespaceMagics


# In[2]:

class VariableInspectorWindow(object):
    instance = None
    
    def __init__(self, ipython):
        """Public constructor."""
        if VariableInspectorWindow.instance is not None:
            raise Exception("""Only one instance of the Variable Inspector can exist at a 
                time.  Call close() on the active instance before creating a new instance.
                If you have lost the handle to the active instance, you can re-obtain it
                via `VariableInspectorWindow.instance`.""")
        
        VariableInspectorWindow.instance = self
        self.closed = False
        self.namespace = NamespaceMagics()
        self.namespace.shell = ipython.kernel.shell
        
        self._popout = widgets.PopupWidget()
        self._popout.description = "Variable Inspector"
        self._popout.button_text = self._popout.description

        self._modal_body = widgets.ContainerWidget()
        self._modal_body.set_css('overflow-y', 'scroll')

        self._modal_body_label = widgets.HTMLWidget(value = 'Not hooked')
        self._modal_body.children = [self._modal_body_label]

        self._popout.children = [
            self._modal_body, 
        ]
        
        self._ipython = ipython
        self._ipython.register_post_execute(self._fill)
        
    def close(self):
        """Close and remove hooks."""
        if not self.closed:
            del self._ipython._post_execute[self._fill]
            self._popout.close()
            self.closed = True
            VariableInspectorWindow.instance = None

    def _fill(self):
        """Fill self with variable information."""
        values = self.namespace.who_ls()
        self._modal_body_label.value = '<table class="table table-bordered table-striped"><tr><th>Name</th><th>Type</th><th>Value</th></tr><tr><td>'+'</td></tr><tr><td>'.join(['{0}</td><td>{1}</td><td>{2}'.format(v, type(eval(v)).__name__, str(eval(v))) for v in values])+'</td></tr></table>'

    def _ipython_display_(self):
        """Called when display() or pyout is used to display the Variable 
        Inspector."""
        self._popout._ipython_display_()
        self._popout.add_class('vbox')
        self._modal_body.add_class('box-flex1')


# In[3]:

inspector = VariableInspectorWindow(get_ipython())
inspector

# In[ ]:



