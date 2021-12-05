##-----------------------------------------------------------------------------
## Primitive: Template
##-----------------------------------------------------------------------------
# class Template(BasePrimitive):
#     """
#     """
#     def __init__(self, action, context):
#         BasePrimitive.__init__(self, action, context)
#         self.log = context.pipeline_logger
#         self.cfg = self.context.config.instrument
# 
#     def _pre_condition(self):
#         """Check for conditions necessary to run this process"""
#         checks = []
#         return np.all(checks)
# 
#     def _post_condition(self):
#         """Check for conditions necessary to verify that the process run correctly"""
#         checks = []
#         return np.all(checks)
# 
#     def _perform(self):
#         """
#         Returns an Argument() with the parameters that depends on this operation.
#         """
#         self.log.info(f"Running {self.__class__.__name__} action")
# 
#         return self.action.args

