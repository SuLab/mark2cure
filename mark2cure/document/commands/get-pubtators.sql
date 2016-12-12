SELECT  `document_pubtator`.`document_id`,
        ANY_VALUE(`document_pubtator`.`content`),
        GROUP_CONCAT(DISTINCT `document_section`.`id`) as `section_ids`

FROM `document_pubtator`

JOIN `document_section`
    ON `document_section`.`document_id` = `document_pubtator`.`document_id`

WHERE `document_pubtator`.`content` != ''
  AND `document_pubtator`.`document_id` IN ({0})

GROUP BY `document_pubtator`.`document_id`;
