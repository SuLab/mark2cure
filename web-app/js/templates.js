define(function(require) {
  'use strict';
  return {
    blocks : {
      header      : require('tpl!templates/blocks/header.tmpl'),
      footer      : require('tpl!templates/blocks/footer.tmpl'),
    },

    modals : {
      first_run   : require('tpl!templates/modals/first_run.tmpl'),
      complete    : require('tpl!templates/modals/complete.tmpl'),

      survey      : require('tpl!templates/modals/survey.tmpl'),
      message     : require('tpl!templates/modals/message.tmpl'),

      settings      : require('tpl!templates/modals/settings.tmpl'),
      instructions  : require('tpl!templates/modals/instructions.tmpl')
    },

    library : {
      index   : require('tpl!templates/library/index.tmpl'),
      item    : require('tpl!templates/library/item.tmpl')
    },

    paragraph : {
      index     : require('tpl!templates/paragraph/index.tmpl'),
      worditem  : require('tpl!templates/paragraph/worditem.tmpl')
    },

    analytics : {
      network : require('tpl!templates/analytics/network.tmpl')
    },

    snippets : {
      paragraph_info  : require('tpl!templates/snippets/paragraph_info.tmpl'),
      network_info    : require('tpl!templates/snippets/network_info.tmpl')
    }

  };
});
