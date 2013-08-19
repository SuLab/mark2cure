define(['backbone', 'models/User'], 
    
        function(Backbone, User) {
  'use strict';

  return Backbone.Collection.extend({
    model   : User,
    url     : '/api/v1/user',
  });
});
