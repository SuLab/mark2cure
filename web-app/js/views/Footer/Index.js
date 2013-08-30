define(['marionette', 'templates'], 

        function (Marionette, templates) {
  'use strict';

  return Marionette.ItemView.extend({
    template : templates.footer.index,
    className : 'container'
  });
});
